#Step 0: Libraries and paths#

library(dplyr)
library(sesame)
library(sesameData)
library(GenomicRanges)
library(sva)
library(ggplot2)
library(BiocParallel)
library(Rtsne)
library(randomcoloR)
register(SnowParam(8))


## Directories (unchanged)
idat_dir_ref   <- "~/MINT_Meth/HeidelbergRef"
idat_dir_lucas <- "~/MINT_Meth/Lucas"
idat_dir_mint  <- "~/MINT_Meth/MINT"

idat_dirs <- c(idat_dir_ref, idat_dir_lucas, idat_dir_mint)

output_dir <- "~/MINT_Meth/RObjects/PreprocessingSesame/QC"

## Samplesheet
samplesheet_path <- "~/MINT_Meth/Samplesheets/SamplesheetCapper.csv"
samplesheet <- read.csv(samplesheet_path, stringsAsFactors = FALSE)
samplesheet$Basename_full <- sub(
  "_Grn\\.idat$",
  "",
  samplesheet$File_grn
)

## Annotation
##manifest_path <- "~/MINT_Meth/RObjects/tSNE/Full_dataset/SeSAMe_all/openseSAMe/EPIC-8v2-0_A1.csv"
##manifest <- read.csv(manifest_path, stringsAsFactors = FALSE)

hm450_anno <- as.data.frame(
  sesameData_getManifestGRanges("HM450")
)

## Probe filter lists (unchanged)
Zhou_input          <- "~/MINT_Meth/Probe_Filter/AppendixD_Zhou_et_al_MASKgeneral_list.txt"
CrossReactive_input <- "~/MINT_Meth/Probe_Filter/AppendixE_CrossReactiveProbes_EPICv1.txt"
Problematic_input   <- "~/MINT_Meth/Probe_Filter/AppendixF_ProblematicProbes_EPICv1-b5.txt"

Zhou         <- read.table(Zhou_input, stringsAsFactors = FALSE)[,1]
CrossReactive <- read.table(CrossReactive_input, stringsAsFactors = FALSE)[,1]
Problematic   <- read.table(Problematic_input, stringsAsFactors = FALSE)[,1]

remove_probes <- unique(c(Zhou, CrossReactive, Problematic))

# ============================
# STEP 1: Find and Match IDATs
# ============================

idat_prefixes <- unlist(
  lapply(
    idat_dirs,
    searchIDATprefixes
  )
)

idat_df <- data.frame(
  prefix = idat_prefixes,
  Filename = basename(idat_prefixes),
  stringsAsFactors = FALSE
) %>%
  left_join(
    samplesheet,
    by = "Filename"
  )

stopifnot(!any(is.na(idat_df$Sample_ID)))

stopifnot(
  length(unique(idat_df$Sample_ID)) ==
    nrow(idat_df)
)

# ============================
# STEP 2: Read IDATs and perform Platform Check
# ============================

sdfs <- bplapply(
  idat_df$prefix,
  readIDATpair
)

names(sdfs) <- idat_df$Sample_ID

cat(
  "Samples loaded:",
  length(sdfs),
  "\n"
)

saveRDS(
  sdfs,
  file.path(output_dir, "sdfs_preQC.rds")
)

probe_counts <- sapply(
  sdfs,
  nrow
)

platform_df <- data.frame(
  Sample_ID = names(probe_counts),
  n_probes = probe_counts
)

write.csv(
  platform_df,
  file.path(output_dir, "Platform_probe_counts.csv"),
  row.names = FALSE
)

##############################################################

# ============================
# STEP 3: Raw QC Metrics
# ============================

qc_list <- bplapply(
  sdfs,
  sesameQC_calcStats
)

names(qc_list) <- names(sdfs)

qc_df <- bind_rows(
  lapply(
    seq_along(qc_list),
    function(i){
      
      data.frame(
        Sample_ID = names(qc_list)[i],
        as.data.frame(as.list(qc_list[[i]]@stat)),
        check.names = FALSE
      )
      
    }
  )
)

qc_df <- qc_df %>%
  left_join(
    samplesheet,
    by = "Sample_ID"
  )

write.csv(
  qc_df,
  file.path(output_dir, "QC_raw_metrics.csv"),
  row.names = FALSE
)

rm(qc_list)
gc()

# ============================
# STEP 4: Intensity QC
# ============================

intensity_df <- data.frame(
  Sample_ID = names(sdfs),
  intensity = sapply(
    sdfs,
    function(x){
      
      mean(
        x$MG +
          x$MR +
          x$UG +
          x$UR,
        na.rm = TRUE
      )
      
    }
  )
)

intensity_df <- intensity_df %>%
  left_join(
    samplesheet,
    by = "Sample_ID"
  )

write.csv(
  intensity_df,
  file.path(output_dir, "Signal_Intensity.csv"),
  row.names = FALSE
)

############################################################

# ============================
# STEP 5: Build QC Summary Table
# ============================
qc_summary <- qc_df %>%
  select(
    Sample_ID,
    frac_dt,
    frac_dt_cg,
    mean_intensity,
    medR,
    medG
  ) %>%
  left_join(
    intensity_df,
    by = "Sample_ID"
  ) %>%
  left_join(
    samplesheet,
    by = "Sample_ID"
  )

write.csv(
  qc_summary,
  file.path(
    output_dir,
    "QC_Summary_AllSamples.csv"
  ),
  row.names = FALSE
)

# ============================
# STEP 6: QCDPB preprocessing (Sesame standard) and beta extraction
# ============================

sdfs_proc <- bplapply(
  sdfs,
  prepSesame,
  prep = "QCDPB"
)

saveRDS(
  sdfs_proc,
  "~/MINT_Meth/RObjects/PreprocessingSesame/sdfs_proc.rds"
)

betas <- bplapply(
  sdfs_proc,
  getBetas,
  collapseToPfx = TRUE
)

rm(sdfs_proc)
gc()

#################################################################

# ============================
# STEP 7: Harmonize betas to 450K using mLiftOver
# ============================

betas <- lapply(
  betas,
  function(x) mLiftOver(x, "HM450")
)

saveRDS(
  betas,
  "~/MINT_Meth/RObjects/PreprocessingSesame/betas_lifted.rds"
)

# ============================
# STEP 8: Intersect CpG probes and create beta matrix 
# ============================

common_probes <- Reduce(
  intersect,
  lapply(betas, names)
)

beta_matrix <- do.call(
  cbind,
  lapply(
    betas,
    function(x) x[common_probes]
  )
)

rownames(beta_matrix) <- common_probes
colnames(beta_matrix) <- names(betas)

rm(betas)
gc()

# ============================
# STEP 9: Probe filtering (recommended + sex chromosomes)
# ============================

keep_probes <- setdiff(
  rownames(beta_matrix),
  remove_probes
)

beta_matrix <- beta_matrix[
  keep_probes,
  ,
  drop = FALSE
]

anno_beta <- hm450_anno[
  match(
    rownames(beta_matrix),
    rownames(hm450_anno)
  ),
]

keep_chr <-
  !is.na(anno_beta$seqnames) &
  !(anno_beta$seqnames %in%
      c("chrX","chrY","X","Y"))

beta_matrix <- beta_matrix[
  keep_chr,
  ,
  drop = FALSE
]

saveRDS(
  beta_matrix,
  "~/MINT_Meth/RObjects/PreprocessingSesame/beta_matrix_filtered.rds"
)

# ============================
# STEP 10: CREATE METADATA
# ============================

metadata_final <- samplesheet[
  match(
    colnames(beta_matrix),
    samplesheet$Sample_ID
  ),
]

saveRDS(
  metadata_final,
  "~/MINT_Meth/RObjects/PreprocessingSesame/metadata_final.rds"
)

# ============================
# STEP X: REMOVE SAMPLES AFTER QC (Not in primary pipeline)
# ============================

# ============================
# STEP 12: Optional Sample Exclusion
# ============================

# bad_samples <- c(
#   "M_16",
#   "M_21",
#   "M_03"
# )

# if(length(bad_samples) > 0){

#   keep_samples <- setdiff(
#     colnames(beta_matrix),
#     bad_samples
#   )

#   beta_matrix <- beta_matrix[
#     ,
#     keep_samples,
#     drop = FALSE
#   ]

#   metadata_final <- metadata_final[
#     metadata_final$Sample_ID %in% keep_samples,
#   ]

#   metadata_final <- metadata_final[
#     match(
#       colnames(beta_matrix),
#       metadata_final$Sample_ID
#     ),
#   ]

#   stopifnot(
#     all(
#       colnames(beta_matrix) ==
#         metadata_final$Sample_ID
#     )
#   )

# }

# cat(
#   "Samples after QC review:",
#   ncol(beta_matrix),
#   "\n"
# )

# ============================
# STEP 11: Dataset Integrity Check
# ============================

sample_na <- colMeans(is.na(beta_matrix))

qc_final <- data.frame(
  Metric = c(
    "Samples",
    "Probes",
    "Total_NA",
    "Median_sample_missingness",
    "Max_sample_missingness"
  ),
  Value = c(
    ncol(beta_matrix),
    nrow(beta_matrix),
    sum(is.na(beta_matrix)),
    median(sample_na),
    max(sample_na)
  )
)

print(qc_final)

cohort_counts <- metadata_final %>%
  count(Cohort)

print(cohort_counts)

write.csv(
  qc_final,
  file.path(
    output_dir,
    "Final_Dataset_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  cohort_counts,
  file.path(
    output_dir,
    "Final_Cohort_Counts.csv"
  ),
  row.names = FALSE
)

write.csv(
  data.frame(
    Sample_ID = names(sample_na),
    frac_na = sample_na
  ),
  file.path(
    output_dir,
    "Sample_Missingness.csv"
  ),
  row.names = FALSE
)

stopifnot(
  length(sdfs) == ncol(beta_matrix)
)

stopifnot(
  all(
    colnames(beta_matrix) ==
      metadata_final$Sample_ID
  )
)

###########################################################
# ============================
# STEP 12: TOP VARIABLE PROBES
# ============================

v <- apply(
  beta_matrix,
  1,
  sd,
  na.rm = TRUE
)

top <- names(
  sort(
    v,
    decreasing = TRUE
  )
)[1:35000]

beta_use_nobc <- beta_matrix[
  top,
  ,
  drop = FALSE
]

saveRDS(
  beta_use_nobc,
  "~/MINT_Meth/RObjects/PreprocessingSesame/beta_use_nobc.rds"
)

# ============================
# STEP 13: PCA
# ============================
pca_nobc <- prcomp(
  t(beta_use_nobc),
  center = TRUE,
  scale. = FALSE
)

saveRDS(
  pca_nobc,
  "~/MINT_Meth/RObjects/PreprocessingSesame/pca_nobc.rds"
)

var_expl <- pca_nobc$sdev^2 /
  sum(pca_nobc$sdev^2)

pc1_lab <- paste0(
  "PC1 (",
  round(var_expl[1] * 100, 1),
  "%)"
)

pc2_lab <- paste0(
  "PC2 (",
  round(var_expl[2] * 100, 1),
  "%)"
)

df <- data.frame(
  Sample_ID = rownames(pca_nobc$x),
  PC1 = pca_nobc$x[,1],
  PC2 = pca_nobc$x[,2]
)

df <- left_join(
  df,
  metadata_final,
  by = "Sample_ID"
)

shape_map <- c(
  "MINT" = 18,
  "Lucas" = 17,
  "Ref" = 16
)

df$size <- ifelse(
  df$Cohort == "MINT",
  4.5,
  3.0
)

df$stroke <- ifelse(
  df$Cohort == "MINT",
  0.6,
  0
)

ggplot(
  df,
  aes(
    PC1,
    PC2,
    color = Methylation_Class,
    shape = Cohort
  )
) +
  geom_point(
    aes(
      size = size,
      stroke = stroke
    ),
    alpha = 0.9
  ) +
  scale_shape_manual(
    values = shape_map
  ) +
  scale_size_identity() +
  xlab(pc1_lab) +
  ylab(pc2_lab) +
  ggtitle("PCA_nobc") +
  theme_minimal(base_size = 14)


########################################################
## Step 14: tSNE ON PCA WITHOUT COMBAT
########################################################

set.seed(1)

tsnepca_nobc <- Rtsne(
  pca_nobc$x[, 1:50],   # use first 50 PCs
  dims = 2,
  perplexity = 28,
  theta = 0,
  max_iter = 5000
)

saveRDS(
  tsnepca_nobc,
  file = "~/MINT_Meth/RObjects/PreprocessingSesame/tsnepca_nobc.rds"
)

########################################################
## Step 15: PLOT tSNE WITHOUT COMBAT
########################################################

df <- data.frame(
  tSNE1 = tsnepca_nobc$Y[,1],
  tSNE2 = tsnepca_nobc$Y[,2],
  Methylation_Class = metadata_final$Methylation_Class,
  Cohort = metadata_final$Cohort
)

## extract coordinates
x <- tsnepca_nobc$Y[,2]
y <- -tsnepca_nobc$Y[,1]

shape_map <- c(
  Ref = 16,
  Lucas = 17,
  MINT = 18
)

pch_vec <- shape_map[metadata_final$Cohort]

lucas_idx <- metadata_final$Cohort == "Lucas"
mint_idx  <- metadata_final$Cohort == "MINT"

class_labels <- metadata_final$Methylation_Class
class_labels[class_labels %in% c("", " ", "NA", "NA ")] <- NA

unique_classes <- sort(
  unique(class_labels),
  na.last = TRUE
)

unique_classes <- ifelse(
  is.na(unique_classes),
  "NA",
  unique_classes
)

## COLORS (UNCHANGED)
color_map <- setNames(
  distinctColorPalette(length(unique_classes)),
  unique_classes
)

color_map["NA"] <- "grey70"

tmp_labels <- class_labels
tmp_labels[is.na(tmp_labels)] <- "NA"

cols_vec <- color_map[tmp_labels]

########################################################
## MAIN PANEL
########################################################

par(mar = c(5, 5, 4, 2))
par(fig = c(0, 0.70, 0, 1))

plot(
  x,
  y,
  col = cols_vec,
  pch = pch_vec,
  main = "tSNEpca_nobc: methylation class (color) + cohort (shape)",
  xlab = "tSNE1",
  ylab = "tSNE2",
  asp = 1,
  xlim = c(-50,20),
  ylim = c(10,60)
)

########################################################
## HIGHLIGHT MINT
########################################################

points(
  x[mint_idx],
  y[mint_idx],
  pch = 18,
  col = cols_vec[mint_idx],
  cex = 1.3
)

points(
  x[mint_idx],
  y[mint_idx],
  pch = 5,
  col = "black",
  cex = 1.4,
  lwd = 0.7
)

########################################################
## HIGHLIGHT LUCAS
########################################################

points(
  x[lucas_idx],
  y[lucas_idx],
  pch = 17,
  col = cols_vec[lucas_idx],
  cex = 1.1
)

points(
  x[lucas_idx],
  y[lucas_idx],
  pch = 2,
  col = "black",
  cex = 1.3,
  lwd = 0.7
)

########################################################
## LABELS (Lucas + MINT)
########################################################

text(
  x[lucas_idx | mint_idx],
  y[lucas_idx | mint_idx],
  labels = metadata_final$Sample_ID[lucas_idx | mint_idx],
  pos = 3,
  cex = 0.6,
  col = "black"
)

########################################################
## LEGENDS
########################################################

par(fig = c(0.70, 1, 0.65, 1), new = TRUE)
par(mar = c(0, 0, 0, 0))
plot.new()

legend(
  "topleft",
  legend = c("Ref", "Lucas", "MINT"),
  pch = c(16, 17, 15),
  col = "black",
  title = "Cohort",
  bty = "n",
  cex = 0.7
)

par(fig = c(0.70, 1, 0.15, 0.75), new = TRUE)
plot.new()

legend(
  "bottomleft",
  legend = unique_classes,
  pch = 16,
  col = color_map[unique_classes],
  title = "Methylation Class",
  bty = "n",
  cex = 0.75,
  ncol = 2,
  y.intersp = 0.5
)






















# ============================
# STEP X1: PCA AND tSNE WITH COMBAT (sensitivity_)
# ============================
batch <- factor(
  metadata_final$Array_Cluster
)

stopifnot(
  length(batch) ==
    ncol(mval_matrix)
)
table(batch)
sum(is.na(mval_matrix))
#############################################

mval_bc <- ComBat(
  dat = mval_matrix,
  batch = batch,
  mod = NULL
)

beta_bc <- minfi::ilogit2(mval_bc)

saveRDS(
  mval_bc,
  file="~/MINT_Meth/RObjects/PreprocessingSesame/mval_bc.rds"
)

saveRDS(
  beta_bc,
  file="~/MINT_Meth/RObjects/PreprocessingSesame/beta_bc.rds"
)

beta_use_bc <- beta_bc[
  top,
  ,
  drop = FALSE
]

saveRDS(
  beta_use_bc,
  "~/MINT_Meth/RObjects/PreprocessingSesame/beta_use_bc.rds"
)

pca_bc <- prcomp(
  t(beta_use_bc),
  center = TRUE,
  scale. = FALSE
)

saveRDS(
  pca_bc,
  "~/MINT_Meth/RObjects/PreprocessingSesame/pca_bc.rds"
)

var_expl <- pca_bc$sdev^2 /
  sum(pca_bc$sdev^2)

pc1_lab <- paste0(
  "PC1 (",
  round(var_expl[1] * 100, 1),
  "%)"
)

pc2_lab <- paste0(
  "PC2 (",
  round(var_expl[2] * 100, 1),
  "%)"
)

df <- data.frame(
  Sample_ID = rownames(pca_bc$x),
  PC1 = pca_bc$x[,1],
  PC2 = pca_bc$x[,2]
)

df <- left_join(
  df,
  metadata_final,
  by = "Sample_ID"
)

shape_map <- c(
  "MINT" = 18,
  "Lucas" = 17,
  "Ref" = 16
)

df$size <- ifelse(
  df$Cohort == "MINT",
  4.5,
  3.0
)

df$stroke <- ifelse(
  df$Cohort == "MINT",
  0.6,
  0
)

ggplot(
  df,
  aes(
    PC1,
    PC2,
    color = Methylation_Class,
    shape = Cohort
  )
) +
  geom_point(
    aes(
      size = size,
      stroke = stroke
    ),
    alpha = 0.9
  ) +
  scale_shape_manual(
    values = shape_map
  ) +
  scale_size_identity() +
  xlab(pc1_lab) +
  ylab(pc2_lab) +
  ggtitle("PCA_bc") +
  theme_minimal(base_size = 14)

#run tSNE with combat
set.seed(1)

tsnepca_bc <- Rtsne(
  pca_bc$x[, 1:50],   # use first 50 PCs
  dims = 2,
  perplexity = 28,
  theta = 0,
  max_iter = 5000
)

saveRDS(
  tsnepca_bc,
  file = "~/MINT_Meth/RObjects/PreprocessingSesame/tsnepca_bc.rds"
)

#Plot tSNE with combat
df <- data.frame(
  tSNE1 = tsnepca_bc$Y[,1],
  tSNE2 = tsnepca_bc$Y[,2],
  Methylation_Class = metadata_final$Methylation_Class,
  Cohort = metadata_final$Cohort
)

## extract coordinates
x <- tsnepca_bc$Y[,1]
y <- tsnepca_bc$Y[,2]

shape_map <- c(
  Ref = 16,
  Lucas = 17,
  MINT = 18
)

pch_vec <- shape_map[metadata_final$Cohort]

lucas_idx <- metadata_final$Cohort == "Lucas"
mint_idx  <- metadata_final$Cohort == "MINT"

class_labels <- metadata_final$Methylation_Class
class_labels[class_labels %in% c("", " ", "NA", "NA ")] <- NA

unique_classes <- sort(
  unique(class_labels),
  na.last = TRUE
)

unique_classes <- ifelse(
  is.na(unique_classes),
  "NA",
  unique_classes
)

## COLORS (UNCHANGED)
color_map <- setNames(
  distinctColorPalette(length(unique_classes)),
  unique_classes
)

color_map["NA"] <- "grey70"

tmp_labels <- class_labels
tmp_labels[is.na(tmp_labels)] <- "NA"

cols_vec <- color_map[tmp_labels]

########################################################
## MAIN PANEL
########################################################

par(mar = c(5, 5, 4, 2))
par(fig = c(0, 0.70, 0, 1))

plot(
  x,
  y,
  col = cols_vec,
  pch = pch_vec,
  main = "tSNEpca_bc",
  xlab = "tSNE1",
  ylab = "tSNE2",
  asp = 1,
  xlim = c(-60,-10),
  ylim = c(-15,20)
)

########################################################
## HIGHLIGHT MINT
########################################################

points(
  x[mint_idx],
  y[mint_idx],
  pch = 18,
  col = cols_vec[mint_idx],
  cex = 1.3
)

points(
  x[mint_idx],
  y[mint_idx],
  pch = 5,
  col = "black",
  cex = 1.4,
  lwd = 0.7
)

########################################################
## HIGHLIGHT LUCAS
########################################################

points(
  x[lucas_idx],
  y[lucas_idx],
  pch = 17,
  col = cols_vec[lucas_idx],
  cex = 1.1
)

points(
  x[lucas_idx],
  y[lucas_idx],
  pch = 2,
  col = "black",
  cex = 1.3,
  lwd = 0.7
)

########################################################
## LABELS (Lucas + MINT)
########################################################

text(
  x[lucas_idx | mint_idx],
  y[lucas_idx | mint_idx],
  labels = metadata_final$Sample_ID[lucas_idx | mint_idx],
  pos = 3,
  cex = 0.6,
  col = "black"
)

########################################################
## LEGENDS
########################################################

par(fig = c(0.70, 1, 0.65, 1), new = TRUE)
par(mar = c(0, 0, 0, 0))
plot.new()

legend(
  "topleft",
  legend = c("Ref", "Lucas", "MINT"),
  pch = c(16, 17, 15),
  col = "black",
  title = "Cohort",
  bty = "n",
  cex = 0.7
)

par(fig = c(0.70, 1, 0.15, 0.75), new = TRUE)
plot.new()

legend(
  "bottomleft",
  legend = unique_classes,
  pch = 16,
  col = color_map[unique_classes],
  title = "Methylation Class",
  bty = "n",
  cex = 0.75,
  ncol = 2,
  y.intersp = 0.5
)

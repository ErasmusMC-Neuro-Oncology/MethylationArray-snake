##Preprocess .idat files and create tSNE for MINT-study (joint normalization)
##Author: Lara Orlandini (l.orlandini@erasmusmc.nl)

#History
# 11-05-2026: File creation

##==========================
## 0. PACKAGES and PATHWAYS TO IDATs, SAMPLE SHEET and RECOMMENDED PROBE FILTERS
##==========================

library(minfi)
library(Rtsne)
library(RColorBrewer)
library(limma)
library(ggplot2)
library(stringr)
library(readr)
library(dplyr)
library(sva)
library(sesame)
library(sesameData)

library(BiocParallel)
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 16))
register(MulticoreParam(workers = ncores, stop.on.error = FALSE))
register(SnowParam(16))

# run only once per system — SAFE
#sesameDataCache()

if(exists("snakemake")){
    idat_dir_ref     <- snakemake@params[['idat_dir_Capper']]
    idat_dir_lucas   <- snakemake@params[['idat_dir_Lucas']] 
    samplesheet_ref <- snakemake@params[['samplesheet_Capper']]
    samplesheet_lucas <- snakemake@params[['samplesheet_Lucas']]

    Zhou_input <- snakemake@params[['Zhou_probes']]
    CrossReactive_input <- snakemake@params[['CrossReactive_probes']]
    Problematic_input <- snakemake@params[['Problematic_probes']]
    samplesheet_mint <-snakemake@input[[1]]
}else{
    idat_dir_ref <- '~/mnt/BIGR_home/MINT/data/idat/Capper'
    idat_dir_lucas <- '~/mnt/BIGR_home/MINT/data/idat/Lucas'
    idat_dir_mint <- '~/mnt/BIGR_home/MINT/data/idat/'
    samplesheet_mint <- '~/mnt/BIGR_home/MINT/output/Methylation/samplesheets/Samplesheet_Methylation.csv'
    samplesheet_ref <- '~/mnt/BIGR_home/MINT/data/samplesheet_capper.csv'
    samplesheet_lucas <- '~/mnt/BIGR_home/MINT/data/SampleSheetLucas.csv'

    samplesheet_tsne_lucas <- '~/mnt/BIGR_home/MINT/data/Samplesheet_tSNELucas.csv'
    
    
    Zhou_input <-  '~/mnt/BIGR_home/Resources/EPIC/manifest/AppendixD_Zhou_et_al_MASKgeneral_list.txt' 
    CrossReactive_input <- '~/mnt/BIGR_home/Resources/EPIC/manifest/AppendixE_CrossReactiveProbes_EPICv1.txt' 
    Problematic_input <-  '~/mnt/BIGR_home/Resources/EPIC/manifest/AppendixF_ProblematicProbes_EPICv1-b5.txt' 
}


samplesheet_MINT <- read.csv(samplesheet_mint, stringsAsFactors = FALSE)
samplesheet_ref <- read.csv(samplesheet_ref)
samplesheet_lucas <- read.csv(samplesheet_lucas)

samplesheet <- rbind(
    samplesheet_MINT %>% mutate(
                        Methylation_Class = NA,
                        Cohort = 'MINT'
                    ),
    samplesheet_lucas %>% mutate(patient = paste0('L_',ID), sample = patient,
                                 idat_green = paste0(idat_dir_lucas, '/',basename(idat_green)), 
                                 idat_red = paste0(idat_dir_lucas, '/',basename(idat_red)),
                                 Cohort = 'Lucas') %>%
    select(patient,sample,idat_green,idat_red,Methylation_Class, Cohort),

    samplesheet_ref %>% mutate(
                            idat_green = paste0(idat_dir_ref, '/',basename(idat_green)), 
                            idat_red = paste0(idat_dir_ref, '/',basename(idat_red)),
                            Methylation_Class = Reference.Group.abbreviation,
                            Cohort = 'Capper', Cohort) %>%
    select(patient,sample,idat_green,idat_red,Methylation_Class, Cohort)

) %>%
    mutate(Basename_full = sub("_Grn\\.idat$|\\.gz", "", idat_green),
           Basename_full = gsub('_Grn\\.idat','',Basename_full),
           Basename_full = gsub('/trinity/home/r115502/','~/',Basename_full)) %>%
    filter(file.exists(paste0(Basename_full, "_Grn.idat"))) %>%
    mutate(Sample_ID = sample)

########################################
### 2. READ DATA
########################################
Zhou <- read.table(Zhou_input, stringsAsFactors = FALSE)[,1]
CrossReactive <- read.table(CrossReactive_input, stringsAsFactors = FALSE)[,1]
Problematic <- read.table(Problematic_input, stringsAsFactors = FALSE)[,1]
remove_probes <- unique(c(Zhou, CrossReactive, Problematic))

hm450_anno <- as.data.frame(
  sesameData_getManifestGRanges("HM450")
)

# ============================
# STEP 2: Read IDATs and perform Platform Check
# ============================
if(!file.exists('sdfs.Rds')){
    sdfs <- bplapply(
        samplesheet$Basename_full ,
        readIDATpair
    )
    saveRDS(sdfs,'sdfs.Rds')
    names(sdfs) <- samplesheet$sample
}else{
    sdfs <- readRDS('sdfs.Rds')
}

probe_counts <- sapply(
  sdfs,
  nrow
)

platform_df <- data.frame(
  Sample_ID = samplesheet$sample,
  n_probes = probe_counts
)



write.csv(
  platform_df,
  file.path("Platform_probe_counts.csv"),
  row.names = FALSE
)

##############################################################

# ============================
# STEP 6: QCDPB preprocessing (Sesame standard) and beta extraction
# ============================
if(!file.exists("sdfs_proc.rds")){
    sdfs_proc <- bplapply(
        sdfs,
        prepSesame,
        prep = "QCDPB"
    )
    saveRDS(
        sdfs_proc,
        "sdfs_proc.rds"
    )
}else{
    sdfs_proc <- readRDS("sdfs_proc.rds")
}


# Remove probes 




betas <- bplapply(
  sdfs_proc,
  getBetas,
  collapseToPfx = TRUE,
)

rm(sdfs_proc)
gc()

#################################################################

# ============================
# STEP 7: Harmonize betas to 450K using mLiftOver
# ============================
if(!file.exists("betas_lifted.rds")){
    betas <- lapply(
        betas,
        function(x) mLiftOver(x, "HM450")
    )
    saveRDS(
        betas,
        "betas_lifted.rds"
    )
}else{
    betas <- readRDS("betas_lifted.rds")
}
# ============================
# STEP 8: Intersect CpG probes and create beta matrix 
# ============================

names(betas) <- samplesheet$Sample_ID


common_probes <- Reduce(
  intersect,
  lapply(betas, names)
)

beta <- do.call(
  cbind,
  lapply(
    betas,
    function(x) x[common_probes]
  )
)

rownames(beta) <- common_probes
colnames(beta) <- names(betas)

rm(betas)
gc()

# ============================
# STEP 9: Probe filtering (recommended + sex chromosomes)
# ============================

keep_probes <- setdiff(
  rownames(beta),
  remove_probes
)

beta <- beta[
  keep_probes,
  ,
  drop = FALSE
]

anno_beta <- hm450_anno[
  match(
    rownames(beta),
    rownames(hm450_anno)
  ),
]

keep_chr <-
  !is.na(anno_beta$seqnames) &
  !(anno_beta$seqnames %in%
      c("chrX","chrY","X","Y"))

beta <- beta[
  keep_chr,
  ,
  drop = FALSE
]


saveRDS(
  beta,
  "beta_filtered.rds"
)


beta <- readRDS('beta_filtered.rds')

# ============================
# STEP 10: CREATE METADATA
# ============================

metadata_final <- samplesheet
saveRDS(
  metadata_final,
  "metadata_final.rds"
)
metadata_final <- readRDS('metadata_final.rds')
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
#     colnames(beta),
#     bad_samples
#   )

#   beta <- beta[
#     ,
#     keep_samples,
#     drop = FALSE
#   ]

#   metadata_final <- metadata_final[
#     metadata_final$Sample_ID %in% keep_samples,
#   ]

#   metadata_final <- metadata_final[
#     match(
#       colnames(beta),
#       metadata_final$Sample_ID
#     ),
#   ]

#   stopifnot(
#     all(
#       colnames(beta) ==
#         metadata_final$Sample_ID
#     )
#   )

# }

# cat(
#   "Samples after QC review:",
#   ncol(beta),
#   "\n"
# )

# ============================
# STEP 11: Dataset Integrity Check
# ============================

sample_na <- colMeans(is.na(beta))

qc_final <- data.frame(
  Metric = c(
    "Samples",
    "Probes",
    "Total_NA",
    "Median_sample_missingness",
    "Max_sample_missingness"
  ),
  Value = c(
    ncol(beta),
    nrow(beta),
    sum(is.na(beta)),
    median(sample_na),
    max(sample_na)
  )
)


cohort_counts <- metadata_final %>%
  count(Cohort)

write.csv(
  qc_final,
  file.path(
    
    "Final_Dataset_Summary.csv"
  ),
  row.names = FALSE
)

write.csv(
  cohort_counts,
  file.path(
    
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
    
    "Sample_Missingness.csv"
  ),
  row.names = FALSE
)

stopifnot(
  length(sdfs) == ncol(beta)
)

stopifnot(
  all(
    colnames(beta) ==
      metadata_final$Sample_ID
  )
)

###########################################################
# ============================
# STEP 12: TOP VARIABLE PROBES
# ============================


v <- apply(
  beta,
  1,
  sd,
  na.rm = TRUE
)

top <- names(
  sort(
    v,
    decreasing = TRUE
  )
)[1:32000]

beta_use_nobc <- beta[
  top,
  ,
  drop = FALSE
]

saveRDS(
  beta_use_nobc,
  "beta_use_nobc.rds"
)

# ============================
# STEP 13: PCA
# ============================
beta_use_nobc <- readRDS("beta_use_nobc.rds")
beta_use_nobc <- beta_use_nobc[complete.cases(beta_use_nobc),]


pca_nobc <- prcomp(
  t(beta_use_nobc),
  center = TRUE,
  scale. = FALSE
)

saveRDS(
  pca_nobc,
  "pca_nobc.rds"
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


metadata_final %>% filter(is.na(Cohort))

df <- data.frame(
  Sample_ID = rownames(pca_nobc$x),
  PC1 = pca_nobc$x[,1],
  PC2 = pca_nobc$x[,2]
)
df <- left_join(
  df,
  metadata_final,
  by = "Sample_ID"
) %>%
    mutate(Cohort = ifelse(is.na(Cohort),'MINT',Cohort))


shape_map <- c(
  "MINT" = 18,
  "Lucas" = 17,
  "Capper" = 16
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

df %>% select(Sample_ID,PC1,PC2, Cohort)

ggplot(
  df,
  aes(
    PC1,
    PC2,
    color = Cohort
    #shape = Cohort
  )
) +
  geom_point(alpha = 0.5) +
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

pca_nobc <- readRDS('pca_nobc.rds')
set.seed(1)
tsnepca_nobc <- Rtsne(
  pca_nobc$x[, 1:50],   # use first 50 PCs
  dims = 2,
  pca=F,
  perplexity = 30,
  theta = 0,
  max_iter = 5000
)

saveRDS(
  tsnepca_nobc,
  file = "tsnepca_nobc_perplexity30_50PCs.rds"
)

########################################################
## Step 15: PLOT tSNE WITHOUT COMBAT
########################################################
library(randomcoloR)
tsnepca_nobc <- readRDS('tsnepca_nobc.rds')
df <- data.frame(
  tSNE1 = tsnepca_nobc$Y[,1],
  tSNE2 = tsnepca_nobc$Y[,2],
  Methylation_Class = metadata_final$Methylation_Class,
  Sample_ID = metadata_final$Sample_ID,
  Cohort = metadata_final$Cohort
)


pdf('tsne_no_bc.pdf', height = 6, width = 7)
df %>% arrange(Cohort) %>% ggplot(aes(tSNE1,tSNE2,color = Cohort)) +
    geom_point() + theme_classic(base_size = 12)
dev.off()

df$Methylation_Class %>% unique()

pdf('tsne_no_bc_class_facet.pdf', height = 6, width = 7*3)
subtype <- df %>% mutate(Is_PA_NF1_associated = Methylation_Class == "PA, NF1-associated",
              PA_class =  case_when(
                  Methylation_Class == "PA, NF1-associated" ~ 'PA, NF1-associated',
                  Methylation_Class == "HGAP" ~ 'HGAP',
                  Methylation_Class == 'ANA_PA' ~ 'Anaplastic PA')) %>%

    ggplot(aes(tSNE1,tSNE2,color = PA_class)) +
    geom_point() + theme_classic(base_size = 12) +
    facet_wrap(~Cohort)
subtype
dev.off()
df$Methylation_Class %>% unique()

pdf('tsne_no_bc_cohort.pdf', height = 6, width = 7)
cohort <- df %>% mutate(Is_PA_NF1_associated = Methylation_Class == "PA, NF1-associated",
              PA_class =  case_when(
                  Methylation_Class == "PA, NF1-associated" ~ 'PA, NF1-associated',
                  Methylation_Class == "HGAP" ~ 'HGAP',
                  Methylation_Class == 'ANA_PA' ~ 'Anaplastic PA')) %>%

    ggplot(aes(tSNE1,tSNE2,color = Cohort)) +
    geom_point() + theme_classic(base_size = 12)
cohort
dev.off()



Classes <- read.csv('/home/jurriaan/mnt/BIGR_home/MINT/data/SamplesheetCapper_with_classes.csv') %>%
    rename(Methylation_class_MINT = Methylation_Class)


library(dplyr)
library(stringr)

df <- df %>%
  mutate(
    m_num = str_extract(Sample_ID, "(?<=_M)\\d+"),
    r_num = str_extract(Sample_ID, "(?<=_R)\\d+"),
    join_key = paste0("M_", m_num, ifelse(r_num == "1", "", paste0("R", r_num)))
  )

result <- df %>%
    left_join(Classes, by = c("join_key" = "Sample_ID")) 





pdf('tsne_no_bc_class_zoom.pdf', height = 6, width = 7)
result %>% filter(tSNE1 > 42, tSNE1 < 80,tSNE2< -5,tSNE2>-30) %>% 
    mutate(Is_PA_NF1_associated = Methylation_Class == "PA, NF1-associated") %>% ggplot(aes(tSNE1,tSNE2,color = Is_PA_NF1_associated)) +
    geom_point() + theme_classic(base_size = 12)
dev.off()

head(df)

Sample_In_clouds <- result %>% filter(tSNE1 > 42, tSNE1 < 80,tSNE2< -5,tSNE2>-30) %>%
    mutate(Class = ifelse(is.na(Methylation_Class) ,Methylation_class_MINT,Methylation_Class))


pca_nobc <- prcomp(
  t(beta_use_nobc[,which(metadata_final$Sample_ID %in% Sample_In_clouds$Sample_ID)]),
  center = TRUE,
  scale. = FALSE
)


data.frame(
    PC1 = pca_nobc$x[,1], PC2 = pca_nobc$x[,2],
    methylation_class = Sample_In_clouds$Class,    
    Cohort = Sample_In_clouds$Cohort.x,
    Sample = Sample_In_clouds$Sample_ID) %>%
    ggplot(aes(PC1, PC2, color = methylation_class, shape = Cohort)) + geom_point() 


plot(pca_nobc$x)

# ============================
# STEP X1: PCA AND tSNE WITH COMBAT (sensitivity_)
# ============================


batch <- factor(
  metadata_final$Cohort
)

#############################################
beta <- readRDS('beta_filtered.rds')
beta <- beta[complete.cases(beta),]

M_matrix <- log2(beta / (1 - beta))

v <- apply(
   M_matrix,
  1,
  sd,
  na.rm = TRUE
)


top <- names(
  sort(
    v,
    decreasing = TRUE
  )
)[1:32000]




dim(M_matrix)

M_matrix[1:10,1:10]
mval_bc <- ComBat(
  dat = M_matrix,
  batch = batch,
  mod = NULL
)
beta_bc <- minfi::ilogit2(mval_bc)




v <- apply(
  beta_bc,
  1,
  sd,
  na.rm = TRUE
)

top <- names(
  sort(
    v,
    decreasing = TRUE
  )
)[1:32000]




saveRDS(
  mval_bc,
  file="mval_bc.rds"
)

saveRDS(
  beta_bc,
  file="beta_bc.rds"
)

beta_use_bc <- beta_bc[
  top,
  ,
  drop = FALSE
]

saveRDS(
  beta_use_bc,
  "beta_use_bc.rds"
)

pca_bc <- prcomp(
  t(beta_use_bc),
  center = TRUE,
  scale. = FALSE
)

saveRDS(
  pca_bc,
  "pca_bc.rds"
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
  file = "tsnepca_bc.rds"
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

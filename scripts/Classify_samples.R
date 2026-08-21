#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Classify_samples.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Classify samples using a pretrained model 
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv:
# Usage: 
#
# TODO:
# 1) 
#
# History:
#  05-05-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
library(minfi)
library(conumee)
library(randomForest)
library(limma)

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input<- snakemake@input[[1]]
    classifier <- snakemake@params[['classifier']]
    CGC_classifier <- snakemake@params[['CGC_classifier']]
    ba_coef <- snakemake@params[['ba_coef']]
    material <- snakemake@params[['material']]
    Rpreprocess <- snakemake@params[['Rpreprocess']]
    filter_dir <- snakemake@params[['filter_dir']]
    CNA_data <- snakemake@params[['CNA_data']]
    output<- snakemake@output[[1]]

}else{
    input <- '/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/Methylation/samplesheets/Samplesheet_Methylation.csv'
    classifier <- '/home/jurriaan/mnt/BIGR_home/MNP_classifier/output/rf.pred.RData'
    CGC_classifier <- '/home/jurriaan/mnt/BIGR_home/Resources/Continuous_Grading_Classifier/assets/CGC-Psi_predictor_probe_based_lm_v1.0_epicv2.Rds'
    ba_coef <- '/home/jurriaan/mnt/BIGR_home/MNP_classifier/output/ba.coef.RData'
    material <- 'FFPE'
    Rpreprocess <- '/home/jurriaan/mnt/BIGR_home/MNP_classifier/scripts/MNPprocessIDAT_functions.R'
    filter_dir <- '/home/jurriaan/mnt/BIGR_home/MNP_classifier/filter/'
    CNA_data <- '/home/jurriaan/mnt/BIGR_home/MNP_classifier/CNV_data/'
    output <- '/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/Methylation/results/Methylation_Classes.txt'
}

# ---------------------------------------------------------------------------
# 1. Read data
# ---------------------------------------------------------------------------
samplesheet <- read.delim(input , sep = ',') 
head(samplesheet)
idat_basename <- gsub('_Grn.idat','',samplesheet$idat_green)


if (!material %in% c("Frozen", "FFPE"))
  stop("Material should be either 'Frozen' or 'FFPE'", call. = FALSE)


source(Rpreprocess)
.isEPICv2_local <- function(object) {
  annotation(object)["array"] == "IlluminaHumanMethylationEPICv2"
}

bgcorrect.illumina.patched <- function(rgSet) {
  if (minfi:::.is450k(rgSet) || minfi:::.isEPIC(rgSet) ||
      .isEPICv2_local(rgSet)) {
    NegControls <- getControlAddress(rgSet, controlType = "NEGATIVE")
  } else if (minfi:::.is27k(rgSet)) {
    NegControls <- getControlAddress(rgSet, controlType = "Negative")
  } else {
    stop("bgcorrect.illumina: unsupported array type")
  }
  Green <- getGreen(rgSet)
  Red   <- getRed(rgSet)
  greenBg <- apply(Green[NegControls, , drop = FALSE], 2,
                   function(x) quantile(x, 0.05))
  redBg   <- apply(Red[NegControls,   , drop = FALSE], 2,
                   function(x) quantile(x, 0.05))
  Green <- pmax(sweep(Green, 2, greenBg, "-"), 1)
  Red   <- pmax(sweep(Red,   2, redBg,   "-"), 1)
  assay(rgSet, "Green") <- Green
  assay(rgSet, "Red")   <- Red
  rgSet
}

# Replace minfi's version in its own namespace so MNPpreprocessIllumina picks it up
assignInNamespace("bgcorrect.illumina",
                  bgcorrect.illumina.patched,
                  ns = "minfi")
# ---------------------------------------------------------------------------
# 2. Referentiebestanden laden
# ---------------------------------------------------------------------------
required <- list(
  file.path(CNA_data, "CNanalysis4_conumee_ANNO.vh20150715.RData"),
  file.path(CNA_data, "CNanalysis4_conumee_REF-M.vh20150715.RData"),
  file.path(CNA_data, "CNanalysis4_conumee_REF-F.vh20150715.RData"),
  file.path(classifier),
  file.path(ba_coef)
)
for (f in required) {
  if (!file.exists(f))
    stop("Ontbrekend bestand: ", f,
         "\nDraai 00_setup.R / 01_preprocessing.R / 02_training.R eerst.",
         call. = FALSE)
}

cat("Laden referentiebestanden...\n")
load(file.path(CNA_data, "CNanalysis4_conumee_ANNO.vh20150715.RData"))  # annoXY
load(file.path(CNA_data, "CNanalysis4_conumee_REF-M.vh20150715.RData")) # refM.data
load(file.path(CNA_data, "CNanalysis4_conumee_REF-F.vh20150715.RData")) # refF.data
load(file.path(classifier))                               # rf.pred
load(file.path(ba_coef))                               # methy.coef, unmethy.coef

# ---------------------------------------------------------------------------
# 3. IDAT inladen en normaliseren
# ---------------------------------------------------------------------------
cat("Inladen IDAT-bestanden...\n")
RGset <- read.metharray(idat_basename, verbose = TRUE, force = TRUE)


if(RGset@annotation[[1]] == 'Unknown'){
    library(IlluminaHumanMethylationEPICv2manifest)    
    annotation(RGset) <- c(
        array = "IlluminaHumanMethylationEPICv2",
        annotation = "ilm10b4.hg19")
    }


cat("MNPpreprocessIllumina normalisatie...\n")
Mset <- MNPpreprocessIllumina(RGset)

# ---------------------------------------------------------------------------
# 5. Probe filtering (zelfde filterlijsten als training)
# ---------------------------------------------------------------------------
cat("\nProbe filtering...\n")
probe_ids <- sub("_.*$", "", rownames(Mset))
amb.filter  <- read.table(file.path(filter_dir, "amb_3965probes.vh20151030.txt"),       header = FALSE)
epic.filter <- read.table(file.path(filter_dir, "epicV1B2_32260probes.vh20160325.txt"), header = FALSE)
snp.filter  <- read.table(file.path(filter_dir, "snp_7998probes.vh20151030.txt"),       header = FALSE)
xy.filter   <- read.table(file.path(filter_dir, "xy_11551probes.vh20151030.txt"),       header = FALSE)
rs.filter   <- grep("rs", rownames(Mset))
ch.filter   <- grep("ch", rownames(Mset))



remove <- unique(c(
  match(amb.filter[, 1],  rownames(Mset)),
  match(epic.filter[, 1], rownames(Mset)),
  match(snp.filter[, 1],  rownames(Mset)),
  match(xy.filter[, 1],   rownames(Mset)),
  rs.filter, ch.filter
))

remove <- remove[!is.na(remove)]
Mset_filtered <- Mset[-remove, ]

# ---------------------------------------------------------------------------
# 6. Batch-effect correctie (conform trainingscoëfficiënten)
# ---------------------------------------------------------------------------
methy   <- getMeth(Mset_filtered)
unmethy <- getUnmeth(Mset_filtered)

# Strip EPICv2 suffixes BEFORE batch correction so coef alignment is correct
# Keep only probes that exist in the training coefficients
rownames(methy)   <- sub("_.*$", "", rownames(methy))
rownames(unmethy) <- sub("_.*$", "", rownames(unmethy))

in_coef <- rownames(methy) %in% names(methy.coef[[material]])
cat("Keeping", sum(in_coef), "of", nrow(methy), "probes for batch correction\n")



methy_bc   <- methy[in_coef, , drop = FALSE]
unmethy_bc <- unmethy[in_coef, , drop = FALSE]

# Align coefficients to row order
coef_methy   <- methy.coef[[material]][rownames(methy_bc)]
coef_unmethy <- unmethy.coef[[material]][rownames(unmethy_bc)]

methy.ba   <- 2^(log2(methy_bc   + 1) + coef_methy)
unmethy.ba <- 2^(log2(unmethy_bc + 1) + coef_unmethy)

betas_sample <- methy.ba / (methy.ba + unmethy.ba + 100)
# rownames are already stripped, skip the sub() that was here before
dim(betas_sample)
# ---------------------------------------------------------------------------
# 7. Tumorclassificatie
# ---------------------------------------------------------------------------

cat("\n--- Tumorclassificatie ---\n")

classifier_probes <- rownames(rf.pred$importance)
available_probes  <- intersect(classifier_probes, rownames(betas_sample))
missing_n         <- length(classifier_probes) - length(available_probes)

if (missing_n > 0)
  warning(missing_n, " classifier-probes ontbreken; geïmputeerd als 0.5", call. = FALSE)

# pred_mat: rows = samples, cols = classifier probes
pred_mat <- matrix(0.5,
                   nrow = ncol(betas_sample),
                   ncol = length(classifier_probes),
                   dimnames = list(colnames(betas_sample), classifier_probes))

# Transpose: betas_sample is probes×samples, pred_mat needs samples×probes
pred_mat[, available_probes] <- t(betas_sample[available_probes, , drop = FALSE])

rf_votes <- predict(rf.pred, newdata = pred_mat, type = "vote")
rf_class  <- colnames(rf_votes)[apply(rf_votes, 1, which.max)]
rf_score  <- apply(rf_votes, 1, max)   # confidence score 0-1

# ---------------------------------------------------------------------------
# 8. Predict CGC
# ---------------------------------------------------------------------------
annotation(RGset)['annotation'] <- "20a1.hg38"

proc <- preprocessNoob(RGset, offset = 0, dyeCorr = T, verbose = TRUE, dyeMethod="single")  #dyeMethod="reference"

mvalue <- ratioConvert(proc, what = "M") |>
  assays() |>
  purrr::pluck('listData') |>
  purrr::pluck("M") |>
  data.table::as.data.table(keep.rownames = "probe_id")
#-------------------------------------------------------------------------------
# 2.2 Predict CGC
#-------------------------------------------------------------------------------
CGC_classifier <- readRDS(CGC_classifier)

# acquire the exact same m-values
data <- mvalue |> 
  tibble::column_to_rownames('probe_id') |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(rownames(CGC_classifier$beta)) |> 
  as.matrix()


# apply lm to the data
CGC <- glmnet::predict.glmnet(CGC_classifier, data)


#-------------------------------------------------------------------------------
# 3.1 Join an write to file
#-------------------------------------------------------------------------------
data.frame(sample = samplesheet$patient,
           class  = rf_class,
           score  = round(rf_score, 3),
           CGC = CGC[,1]) %>%
    write.table(
        output, sep = '\t', row.names = F, quote = F)

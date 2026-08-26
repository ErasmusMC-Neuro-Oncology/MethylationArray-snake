#!/usr/bin/env Rscript
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Preprocess_methylation.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# MULTI-PLATFORM METHYLATION HARMONIZATION PIPELINE
# sesame + mLiftOver | 450K / EPICv1 / EPICv2 -> common HM450 probe space,
# followed by nested ComBat batch correction, written out as a single .h5ad.
#
# Part 1 of 2. Sections 0-12 of the original Methylation_analysis_new.R,
# internals unchanged. Part 2 (Embeddings_Methylation.R) consumes the h5ad.
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv: envs/methylation.yaml
# Usage: snakemake script directive
#
# TODO:
# 1)
#
# History:
#  21-08-2026: Split out of Methylation_analysis_new.R; h5ad output added
#  21-08-2026: Per-sample QC / detection-failure fractions carried into obs
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
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
library(data.table)
library(methylclock)
library(anndata)
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 16))
## SnowParam, not MulticoreParam: preprocessCore's pthread-based quantile
## normalisation fails with "return code from pthread_create() is 22" when
## called inside a forked worker.
register(SnowParam(workers = ncores, stop.on.error = FALSE), default = TRUE)
## reticulate >= 1.41 otherwise provisions its own ephemeral uv-managed Python
## env, ignoring the conda env and its anndata/pandas pins.
reticulate::use_condaenv(Sys.getenv("CONDA_PREFIX"), required = TRUE)
#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if (exists("snakemake")) {
  samplesheet_mint          <- snakemake@input[['samplesheet_mint']]
  idat_dir_ref              <- snakemake@params[['idat_dir_Capper']]
  idat_dir_lucas            <- snakemake@params[['idat_dir_Lucas']]
  samplesheet_ref           <- snakemake@params[['samplesheet_Capper']]
  samplesheet_ref_pediatric <- snakemake@params[['samplesheet_Sturm']]
  samplesheet_lucas         <- snakemake@params[['samplesheet_Lucas']]
  Zhou_input                <- snakemake@params[['Zhou_probes']]
  CrossReactive_input       <- snakemake@params[['CrossReactive_probes']]
  Problematic_input         <- snakemake@params[['Problematic_probes']]
  classes_csv               <- snakemake@params[['classes_csv']]
  OUT_DIR                   <- snakemake@params[['out_dir']]
  adata_out                 <- snakemake@output[['adata']]
} else {
  samplesheet_mint          <- '~/mnt/BIGR_home/MINT/output/Methylation/samplesheets/Samplesheet_Methylation.csv'
  idat_dir_ref              <- '~/mnt/BIGR_home/MINT/data/idat/Capper'
  idat_dir_lucas            <- '~/mnt/BIGR_home/MINT/data/idat/Lucas'
  samplesheet_ref           <- '~/mnt/BIGR_home/MINT/data/samplesheet_capper.csv'
  samplesheet_ref_pediatric <- '~/mnt/BIGR_home/MINT/data/Samplesheet_Sturm.csv'
  samplesheet_lucas         <- '~/mnt/BIGR_home/MINT/data/SampleSheetLucas.csv'
  Zhou_input                <- '~/mnt/BIGR_home/Resources/EPIC/manifest/AppendixD_Zhou_et_al_MASKgeneral_list.txt'
  CrossReactive_input       <- '~/mnt/BIGR_home/Resources/EPIC/manifest/AppendixE_CrossReactiveProbes_EPICv1.txt'
  Problematic_input         <- '~/mnt/BIGR_home/Resources/EPIC/manifest/AppendixF_ProblematicProbes_EPICv1-b5.txt'
  classes_csv               <- '~/mnt/BIGR_home/MINT/data/SamplesheetCapper_with_classes.csv'
  OUT_DIR                   <- '~/mnt/BIGR_home/MINT/output/Methylation/harmonized'
  adata_out                 <- '~/mnt/BIGR_home/MINT/output/Methylation/harmonized/methylation_harmonized.h5ad'
}
# NOTE: samplesheet_ref_pediatric is now a proper parameter. In the original it
# was only defined in the interactive branch, so a snakemake run would have
# failed at Section 1. idat_dir_mint was dropped - it was assigned but never
# read (MINT paths come from the samplesheet's own idat_green column).
INTERMEDIATE_DIR <- file.path(OUT_DIR, "intermediate_chunks")
INTERMEDIATE_DETECTFAIL_DIR <- file.path(OUT_DIR, "intermediate_chunks_detectfail")
LOG_DIR          <- file.path(OUT_DIR, "logs")
dir.create(INTERMEDIATE_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(INTERMEDIATE_DETECTFAIL_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)

## ---- Tunable parameters -----------------------------------------------
CHUNK_SIZE            <- 500
DETECTION_PVAL_THRESH <- 0.01
MAX_SAMPLE_MISSING    <- 0.10
MAX_SAMPLE_FAILED_PROBES <- Inf
MAX_MISSING_FRACTION  <- 0
TARGET_PLATFORM       <- "HM450"
DROP_SEX_CHROM        <- TRUE
WRITE_CSV_EXPORTS     <- FALSE

##==========================
## 1. SAMPLE SHEET
##==========================
samplesheet_MINT  <- read.csv(samplesheet_mint, stringsAsFactors = FALSE)
samplesheet_ref   <- read.csv(samplesheet_ref)
samplesheet_ref_pediatric   <- read.csv(samplesheet_ref_pediatric)
samplesheet_lucas <- read.csv(samplesheet_lucas)

glioma_pattern_ped <- "glioblastoma|astrocytoma|oligodendroglioma|ganglioglioma|glioma|^PA[,/ ]|dysembryoplastic"

selected_classes <- samplesheet_ref_pediatric$meth_class_rf %>%
  unique() %>%
  .[
    (grepl(glioma_pattern_ped, ., ignore.case = TRUE) & 
     !grepl("chordoid glioma", ., ignore.case = TRUE)) |
    grepl("non-classifiable", ., ignore.case = TRUE) |
    grepl("^control", ., ignore.case = TRUE)
  ]
samplesheet_ref %>% head()



samplesheet <- rbind(
    samplesheet_MINT %>% mutate(Methylation_Class = NA, Cohort = 'MINT'),
    
  samplesheet_lucas %>% mutate(
    patient = paste0('L_', ID), sample = patient,
    idat_green = paste0(idat_dir_lucas, '/', basename(idat_green)),
    idat_red   = paste0(idat_dir_lucas, '/', basename(idat_red)),
    Cohort = 'Lucas'
    ) %>% select(patient, sample, idat_green, idat_red, Methylation_Class, Cohort),
  
  samplesheet_ref %>%
  filter(grepl("glioblastoma|astrocytoma|oligodendroglioma|ganglioglioma|glioma|control",Reference.Group..Classifier.version.V11.),
         !grepl("chordoid glioma", Reference.Group..Classifier.version.V11.)) %>%
  mutate(
    idat_green = paste0(idat_dir_ref, '/', basename(idat_green)),
    idat_red   = paste0(idat_dir_ref, '/', basename(idat_red)),
    Methylation_Class = Reference.Group.abbreviation,
    Cohort = 'Capper'
  ) %>% select(patient, sample, idat_green, idat_red, Methylation_Class, Cohort),
  

  samplesheet_ref_pediatric %>% filter(meth_class_rf %in% selected_classes) %>%
  mutate(Methylation_Class = meth_class_rf, patient = Sample_Name, sample = Sample_Name,Cohort = 'Sturm') %>% select(patient,sample,idat_green,idat_red,Methylation_Class,Cohort)

) %>%
  mutate(
    Basename_full = sub("_Grn\\.idat$|\\.gz", "", idat_green),
    Basename_full = gsub('_Grn\\.idat', '', Basename_full),
    Basename_full = gsub('/trinity/home/r115502/', '~/mnt/BIGR_home/', Basename_full)
  ) %>%
  filter(file.exists(paste0(Basename_full, "_Grn.idat"))) %>%
  mutate(Sample_ID = sample)



stopifnot(!any(duplicated(samplesheet$Sample_ID)))
cat(sprintf("[%s] Sample sheet loaded: %d samples\n", Sys.time(), nrow(samplesheet)))

##==========================
## 2. PROBE BLACKLISTS + REFERENCE MANIFEST
##==========================
Zhou           <- read.table(Zhou_input, stringsAsFactors = FALSE)[, 1]
CrossReactive  <- read.table(CrossReactive_input, stringsAsFactors = FALSE)[, 1]
Problematic    <- read.table(Problematic_input, stringsAsFactors = FALSE)[, 1]
remove_probes  <- unique(c(Zhou, CrossReactive, Problematic))
cat(sprintf("[%s] Blacklist probes (Zhou + cross-reactive + problematic): %d\n",
            Sys.time(), length(remove_probes)))

hm450_anno <- as.data.frame(sesameData_getManifestGRanges("HM450"))
hm450_probe_ids <- rownames(hm450_anno)
sexchrom_probes <- rownames(hm450_anno)[hm450_anno$seqnames %in% c("chrX", "chrY")]
cat(sprintf("[%s] HM450 sex-chromosome probes: %d (DROP_SEX_CHROM=%s)\n",
            Sys.time(), length(sexchrom_probes), DROP_SEX_CHROM))

##==========================
## 3. PLATFORM DETECTION
##==========================
detect_platform <- function(prefix) {
  sdf <- tryCatch(sesame::readIDATpair(prefix), error = function(e) NULL)
  if (is.null(sdf)) return(NA_character_)
  plat <- attr(sdf, "platform")
  rm(sdf); gc(FALSE)
  if (is.null(plat) || length(plat) == 0) return(NA_character_)
  as.character(plat)
}

if (!"Platform" %in% colnames(samplesheet)) {
  cat(sprintf("[%s] Reading sesame-detected platform for %d samples (parallel)...\n",
              Sys.time(), nrow(samplesheet)))
  samplesheet$Platform <- unlist(bplapply(
    samplesheet$Basename_full, detect_platform, BPPARAM = bpparam()
  ))

  n_undetected <- sum(is.na(samplesheet$Platform))
  if (n_undetected > 0) {
    warning(sprintf(
      paste0(
        "%d/%d samples had no platform attribute from sesame - falling back ",
        "to probe-count classification for those only."
      ),
      n_undetected, nrow(samplesheet)
    ))
    fallback_by_count <- function(prefix) {
      sdf <- tryCatch(sesame::readIDATpair(prefix), error = function(e) NULL)
      if (is.null(sdf)) return(NA_character_)
      n <- nrow(sdf)
      plat <- if (n < 600000) "HM450" else if (n < 900000) "EPIC" else "EPICv2"
      rm(sdf); gc(FALSE)
      plat
    }
    idx_na <- which(is.na(samplesheet$Platform))
    samplesheet$Platform[idx_na] <- unlist(bplapply(
      samplesheet$Basename_full[idx_na], fallback_by_count, BPPARAM = bpparam()
    ))
  }

  fwrite(samplesheet, file.path(LOG_DIR, "samplesheet_with_platform.csv"))
} else {
  cat(sprintf("[%s] Using existing Platform column from sample sheet.\n", Sys.time()))
}

cat("[Platform distribution]\n")
print(table(samplesheet$Platform, useNA = "ifany"))
samplesheet <- samplesheet %>% filter(Platform %in% c("HM450", "EPIC", "EPICv2"))

##==========================
## 4. PASS 1: per-platform, chunked QC + mLiftOver -> save to disk
##==========================
chunk_registry <- data.frame(
  file = character(),
  platform = character(), n_samples = integer(),
  stringsAsFactors = FALSE
)
sample_qc <- data.frame(
  Sample_ID = character(), Platform = character(),
  n_probes_native = integer(),
  n_failed_native = integer(), frac_missing_native = double(),
  n_failed_hm450 = integer(), frac_missing_hm450 = double(),
  stringsAsFactors = FALSE
)

process_chunk <- function(prefixes, sample_ids, platform, chunk_id) {
  cat(sprintf("[%s]   %s chunk %d: %d samples - openSesame QCDPB (mask=FALSE)...\n",
              Sys.time(), platform, chunk_id, length(prefixes)))

  betas_native <- openSesame(
    prefixes, prep = "QCDPB", func = function(sdf) sesame::getBetas(sdf, mask = FALSE),
    BPPARAM = bpparam()
  )
  if (is.null(dim(betas_native))) {
    betas_native <- matrix(
      betas_native, ncol = 1,
      dimnames = list(names(betas_native), sample_ids)
    )
  } else {
    colnames(betas_native) <- sample_ids
  }

  n_failed_native     <- colSums(is.na(betas_native))
  frac_missing_native <- n_failed_native / nrow(betas_native)

  if (identical(platform, TARGET_PLATFORM)) {
    betas_hm450 <- betas_native[match(hm450_probe_ids, rownames(betas_native)), , drop = FALSE]
    rownames(betas_hm450) <- hm450_probe_ids
  } else {
    betas_hm450_raw <- mLiftOver(
      betas_native, TARGET_PLATFORM, impute = FALSE, BPPARAM = bpparam()
    )
    betas_hm450 <- betas_hm450_raw[match(hm450_probe_ids, rownames(betas_hm450_raw)), , drop = FALSE]
    rownames(betas_hm450) <- hm450_probe_ids
    rm(betas_hm450_raw)
  }
  stopifnot(identical(rownames(betas_hm450), hm450_probe_ids))

  n_failed_hm450     <- colSums(is.na(betas_hm450))
  frac_missing_hm450 <- n_failed_hm450 / nrow(betas_hm450)

  n_probes_native <- nrow(betas_native)

  out_file <- file.path(
    INTERMEDIATE_DIR, sprintf("%s_chunk%03d.rds", platform, chunk_id)
  )
  saveRDS(betas_hm450, out_file, compress = TRUE)

  rm(betas_native, betas_hm450); gc(FALSE)

  list(
    file = out_file,
    platform = platform, n_samples = length(prefixes),
    qc = data.frame(
      Sample_ID = sample_ids, Platform = platform,
      n_probes_native = n_probes_native,
      n_failed_native = n_failed_native,
      frac_missing_native = frac_missing_native,
      n_failed_hm450 = n_failed_hm450,
      frac_missing_hm450 = frac_missing_hm450,
      row.names = NULL
    )
  )
}

required_globals <- c(
  "OUT_DIR", "INTERMEDIATE_DIR", "LOG_DIR",
  "hm450_probe_ids", "remove_probes", "sexchrom_probes",
  "TARGET_PLATFORM", "DETECTION_PVAL_THRESH", "CHUNK_SIZE"
)
missing_globals <- required_globals[!vapply(required_globals, exists, logical(1))]
if (length(missing_globals) > 0) {
  stop("Missing required object(s) before starting the main loop: ",
       paste(missing_globals, collapse = ", "))
}

for (plat in c("HM450", "EPIC", "EPICv2")) {
  plat_sheet <- samplesheet %>% filter(Platform == plat)
  if (nrow(plat_sheet) == 0) next
  n_chunks <- ceiling(nrow(plat_sheet) / CHUNK_SIZE)
  cat(sprintf("[%s] Platform %s: %d samples in %d chunk(s) of <=%d\n",
              Sys.time(), plat, nrow(plat_sheet), n_chunks, CHUNK_SIZE))

  for (i in seq_len(n_chunks)) {
    idx <- ((i - 1) * CHUNK_SIZE + 1):min(i * CHUNK_SIZE, nrow(plat_sheet))
    res <- process_chunk(
      prefixes   = plat_sheet$Basename_full[idx],
      sample_ids = plat_sheet$Sample_ID[idx],
      platform   = plat,
      chunk_id   = i
    )
    chunk_registry <- rbind(chunk_registry, data.frame(
      file = res$file,
      platform = res$platform, n_samples = res$n_samples
    ))
    sample_qc <- rbind(sample_qc, res$qc)
    rm(res); gc(FALSE)
  }
}

fwrite(chunk_registry, file.path(LOG_DIR, "chunk_registry.csv"))
fwrite(sample_qc,      file.path(LOG_DIR, "sample_qc_summary.csv"))
cat(sprintf("[%s] Pass 1 complete: %d chunks written to %s\n",
            Sys.time(), nrow(chunk_registry), INTERMEDIATE_DIR))

##==========================
## 5. SAMPLE-LEVEL QC FILTER
##==========================
fails_fraction_criterion <- sample_qc$frac_missing_native > MAX_SAMPLE_MISSING
fails_count_criterion    <- sample_qc$n_failed_native > MAX_SAMPLE_FAILED_PROBES
bad_samples <- sample_qc$Sample_ID[fails_fraction_criterion | fails_count_criterion]

cat(sprintf(
  "[%s] Excluding %d/%d samples (any-NA fraction>%.0f%%: %d, any-NA count>%s: %d)\n",
  Sys.time(), length(bad_samples), nrow(sample_qc),
  100 * MAX_SAMPLE_MISSING, sum(fails_fraction_criterion),
  ifelse(is.infinite(MAX_SAMPLE_FAILED_PROBES), "Inf (disabled)",
         as.character(MAX_SAMPLE_FAILED_PROBES)),
  sum(fails_count_criterion)
))
kept_samples <- setdiff(sample_qc$Sample_ID, bad_samples)
excluded_qc <- sample_qc[sample_qc$Sample_ID %in% bad_samples, ]
idx_match <- match(excluded_qc$Sample_ID, sample_qc$Sample_ID)
reasons <- Map(function(f, c) {
  hit <- c("any_NA_fraction"[f], "any_NA_count"[c])
  paste(hit[!is.na(hit)], collapse = "+")
}, fails_fraction_criterion[idx_match], fails_count_criterion[idx_match])
excluded_qc$reason <- unlist(reasons)
fwrite(excluded_qc, file.path(LOG_DIR, "excluded_samples_QC.csv"))

##==========================
## 5b. DETECTION-P-VALUE SCAN (per-probe failure count, kept samples only)
##==========================
## "Failed" here means pOOBAH detection p-value > DETECTION_PVAL_THRESH -
## NOT the mask=FALSE NA pattern used elsewhere in this pipeline (that only
## reflects genuinely-unmeasurable signal + mLiftOver structural gaps).
## Raw p-values were never captured during Pass 1 (mask=FALSE via the full
## openSesame() wrapper doesn't expose them), so this re-reads every KEPT
## sample's IDAT pair - a real, additional I/O cost. Kept as cheap as
## possible: only Q/C/D (skips noob entirely, since only p-values are
## needed here, not final betas), and only for samples that survived
## Section 5's QC filter (no point scanning samples we've already dropped).
##
## Uses openSesame()'s own official Q->C->D chain (no manual step-ordering
## risk), then calls pOOBAH() explicitly on that exact processed state -
## reproducing precisely what sesame's internal "P" step would compute if
## it were included in the prep string.
compute_detection_failcount_chunk <- function(prefixes, sample_ids, platform) {
  pvals_native <- openSesame(
    prefixes, prep = "QCD",
    func = function(sdf) sesame::pOOBAH(sdf, return.pval = TRUE),
    BPPARAM = bpparam()
  )
  if (is.null(dim(pvals_native))) {
    pvals_native <- matrix(pvals_native, ncol = 1,
                            dimnames = list(names(pvals_native), sample_ids))
  } else {
    colnames(pvals_native) <- sample_ids
  }

  # sesame convention: a probe FAILS detection when p > threshold (large p
  # = signal indistinguishable from background). NA p-values also count as
  # failed, since detection genuinely couldn't be assessed.
  failed_native <- (pvals_native > DETECTION_PVAL_THRESH) | is.na(pvals_native)
  failed_native <- matrix(as.numeric(failed_native), nrow = nrow(failed_native),
                           dimnames = dimnames(failed_native))

  if (identical(platform, TARGET_PLATFORM)) {
    failed_hm450 <- failed_native[match(hm450_probe_ids, rownames(failed_native)), , drop = FALSE]
  } else {
    # Same row-alignment fix as Pass 1: mLiftOver doesn't guarantee a row
    # for every HM450 target probe, so explicitly reindex its output.
    failed_hm450_raw <- mLiftOver(
      failed_native, TARGET_PLATFORM, impute = FALSE, BPPARAM = bpparam()
    )
    failed_hm450 <- failed_hm450_raw[match(hm450_probe_ids, rownames(failed_hm450_raw)), , drop = FALSE]
    rm(failed_hm450_raw)
  }
  rownames(failed_hm450) <- hm450_probe_ids
  stopifnot(identical(rownames(failed_hm450), hm450_probe_ids))

  # Per-probe failure count contribution from this chunk (a target HM450
  # probe with replicate source probes may get a fractional value here from
  # mLiftOver's averaging - summed as-is; na.rm since structurally-absent
  # probes are NA, not a detection failure, and shouldn't count).
  ## Two summaries from the same scan:
  ##  - per_probe : failure count per HM450 probe, summed over samples (drives
  ##                the Section 6b blacklist extension)
  ##  - per_sample: fraction of that sample's NATIVE probes failing detection.
  ##                Native, not HM450, so it is a property of the array itself
  ##                rather than of the liftover.
  list(
    per_probe  = rowSums(failed_hm450, na.rm = TRUE),
    per_sample = colSums(failed_native, na.rm = TRUE) / nrow(failed_native)
  )
}

n_failed_detection_per_probe <- integer(length(hm450_probe_ids))
names(n_failed_detection_per_probe) <- hm450_probe_ids
frac_failed_detection_per_sample <- numeric(0)

for (plat in c("HM450", "EPIC", "EPICv2")) {
  plat_sheet <- samplesheet %>% filter(Platform == plat, Sample_ID %in% kept_samples)
  if (nrow(plat_sheet) == 0) next
  n_chunks <- ceiling(nrow(plat_sheet) / CHUNK_SIZE)
  cat(sprintf("[%s] Detection-p-value scan %s: %d samples in %d chunk(s)\n",
              Sys.time(), plat, nrow(plat_sheet), n_chunks))
  for (i in seq_len(n_chunks)) {
    idx <- ((i - 1) * CHUNK_SIZE + 1):min(i * CHUNK_SIZE, nrow(plat_sheet))
    chunk_res <- compute_detection_failcount_chunk(
      prefixes   = plat_sheet$Basename_full[idx],
      sample_ids = plat_sheet$Sample_ID[idx],
      platform   = plat
    )
    n_failed_detection_per_probe <- n_failed_detection_per_probe + chunk_res$per_probe
    frac_failed_detection_per_sample <- c(frac_failed_detection_per_sample,
                                          chunk_res$per_sample)
    rm(chunk_res); gc(FALSE)
  }
}
frac_failed_detection_per_probe <- n_failed_detection_per_probe / length(kept_samples)

sample_detection_qc <- data.frame(
  Sample_ID = names(frac_failed_detection_per_sample),
  frac_failed_detection = as.numeric(frac_failed_detection_per_sample),
  row.names = NULL
)
stopifnot(setequal(sample_detection_qc$Sample_ID, kept_samples))
fwrite(sample_detection_qc, file.path(LOG_DIR, "sample_detection_failure.csv"))
cat(sprintf("[%s] Detection-p-value scan complete (%d kept samples)\n",
            Sys.time(), length(kept_samples)))

##==========================
## 6. PASS 2: accumulate per-probe missingness across KEPT samples only
##==========================
n_probes  <- length(hm450_probe_ids)
n_missing_per_probe <- integer(n_probes)
names(n_missing_per_probe) <- hm450_probe_ids
n_kept_total <- 0L

for (i in seq_len(nrow(chunk_registry))) {
  m <- readRDS(chunk_registry$file[i])
  keep_cols <- intersect(colnames(m), kept_samples)
  if (length(keep_cols) == 0) { rm(m); gc(FALSE); next }
  m <- m[, keep_cols, drop = FALSE]
  n_missing_per_probe <- n_missing_per_probe + rowSums(is.na(m))
  n_kept_total <- n_kept_total + ncol(m)
  rm(m); gc(FALSE)
}
stopifnot(n_kept_total == length(kept_samples))

missing_fraction_per_probe <- n_missing_per_probe / n_kept_total


probe_missingness <- data.frame(
  Probe_ID = hm450_probe_ids,
  n_missing = as.integer(n_missing_per_probe),
  n_total = n_kept_total,
  frac_missing = missing_fraction_per_probe,
  n_failed_detection = as.integer(round(n_failed_detection_per_probe)),
  frac_failed_detection = frac_failed_detection_per_probe,
  on_blacklist = hm450_probe_ids %in% remove_probes,
  sex_chrom = hm450_probe_ids %in% sexchrom_probes,
  row.names = NULL
)
fwrite(probe_missingness, file.path(LOG_DIR, "probe_missingness_full.csv"))
cat(sprintf(
  "[%s] Full per-probe missingness table saved (%d probes) -> logs/probe_missingness_full.csv\n",
  Sys.time(), nrow(probe_missingness)
))
print(table(cut(probe_missingness$n_missing,
                 breaks = c(-1, 0, 5, 20, Inf),
                 labels = c("0", "1-5", "6-20", ">20"))))



##==========================
## 6b. EXTEND PROBE BLACKLIST based on DETECTION-P-VALUE failure rate
##==========================
## "Failed" = pOOBAH detection p-value > DETECTION_PVAL_THRESH (from the
## Section 5b scan), NOT the generic mask=FALSE NA pattern (which only
## reflects genuinely-unmeasurable signal + mLiftOver structural gaps and
## says nothing about detection quality). Beyond the curated Zhou/cross-
## reactive/problematic blacklist loaded in Section 2, permanently exclude
## any probe that failed DETECTION in more than this fraction of KEPT
## samples. Deliberately DECOUPLED from MAX_MISSING_FRACTION: if you relax
## MAX_MISSING_FRACTION later to tolerate some missingness in the final
## matrix, badly-behaved probes stay excluded here regardless, rather than
## being let back in under a loosened complete-case tolerance. Reassigns
## `remove_probes` in place, so every downstream use of it (retention
## report, final probe set) picks up the extension automatically with no
## other code changes needed.
PROBE_BLACKLIST_FAIL_FRACTION <- 0.1   # blacklist any probe that FAILED
                                         # DETECTION in more than 10% of
                                         # kept samples


newly_blacklisted <- probe_missingness$Probe_ID[
  probe_missingness$frac_failed_detection > PROBE_BLACKLIST_FAIL_FRACTION
]
n_new_blacklisted <- length(setdiff(newly_blacklisted, remove_probes))
remove_probes <- union(remove_probes, newly_blacklisted)

cat(sprintf(
  paste0(
    "[%s] Extended blacklist: %d probes newly added (FAILED DETECTION in ",
    ">%.1f%% of %d kept samples). Total blacklist now: %d probes.\n"
  ),
  Sys.time(), n_new_blacklisted, 100 * PROBE_BLACKLIST_FAIL_FRACTION,
  length(kept_samples), length(remove_probes)
))
writeLines(remove_probes, file.path(LOG_DIR, "remove_probes_extended.txt"))
cat("  - logs/remove_probes_extended.txt (full blacklist: curated lists + detection-failure additions)\n")

##==========================
## 7. RETENTION REPORT
##==========================
report_thresholds <- c(0, 0.001, 0.01, 0.05)
retention_report <- sapply(report_thresholds, function(th) {
  keep <- names(which(missing_fraction_per_probe <= th))
  keep <- setdiff(keep, remove_probes)
  if (DROP_SEX_CHROM) keep <- setdiff(keep, sexchrom_probes)
  length(keep)
})
retention_df <- data.frame(
  max_missing_fraction = report_thresholds,
  n_probes_retained    = retention_report,
  pct_of_hm450         = round(100 * retention_report / n_probes, 1)
)
cat("\n[Probe retention at different missingness tolerances]\n")
print(retention_df)
fwrite(retention_df, file.path(LOG_DIR, "probe_retention_report.csv"))

##==========================
## 8. FINAL PROBE SET
##==========================
complete_probes <- names(which(missing_fraction_per_probe <= MAX_MISSING_FRACTION))
complete_probes <- setdiff(complete_probes, remove_probes)
if (DROP_SEX_CHROM) complete_probes <- setdiff(complete_probes, sexchrom_probes)
complete_probes <- sort(complete_probes)

cat(sprintf(
  "[%s] Final probe set: %d probes (tolerance=%.3f, after blacklist%s)\n",
  Sys.time(), length(complete_probes), MAX_MISSING_FRACTION,
  if (DROP_SEX_CHROM) " + sex-chrom removal" else ""
))
if (length(complete_probes) == 0) {
  stop("No probes survive filtering at this tolerance.")
}
saveRDS(complete_probes, file.path(LOG_DIR, "final_complete_probes.rds"))

##==========================
## 9. PASS 3: assemble final matrix
##==========================
n_final_samples <- length(kept_samples)
final_matrix <- matrix(
  NA_real_, nrow = length(complete_probes), ncol = n_final_samples,
  dimnames = list(complete_probes, kept_samples)
)


for (i in seq_len(nrow(chunk_registry))) {
  m <- readRDS(chunk_registry$file[i])
  keep_cols <- intersect(colnames(m), kept_samples)
  if (length(keep_cols) == 0) { rm(m); gc(FALSE); next }
  m <- m[complete_probes, keep_cols, drop = FALSE]
  final_matrix[, keep_cols] <- m
  rm(m); gc(FALSE)
}

stopifnot(!anyNA(final_matrix))
cat(sprintf(
  "[%s] Final beta matrix: %d probes x %d samples, 0 missing values (verified)\n",
  Sys.time(), nrow(final_matrix), ncol(final_matrix)
))

##==========================
## 10. SAVE OUTPUTS
##==========================
saveRDS(final_matrix, file.path(OUT_DIR, "beta_matrix_complete_HM450space.rds"))
final_matrix <- readRDS(file.path(OUT_DIR, "beta_matrix_complete_HM450space.rds"))
if (WRITE_CSV_EXPORTS) {
  fwrite(
    data.table(Probe_ID = rownames(final_matrix), final_matrix),
    file.path(OUT_DIR, "beta_matrix_complete_HM450space.csv.gz")
  )
}
fwrite(samplesheet %>% filter(Sample_ID %in% kept_samples),
       file.path(OUT_DIR, "samplesheet_final_kept_samples.csv"))

cat(sprintf("[%s] DONE. Outputs written to %s\n", Sys.time(), OUT_DIR))

##==========================
## 11. BATCH CORRECTION: nested strategy (3-way Step 2)
##     Step 1: WITHIN MINT, batch = Platform (EPIC vs EPICv2)
##     Step 2: Lucas + platform-corrected MINT + Sturm-NF1 subset,
##              batch = Cohort (3 levels)
##     Uncorrected: Capper (all) + Sturm non-NF1 - biologically distinct,
##     carried through unmodified on whatever probe set survives Step 2.
##
## Sturm splits by constitutional NF1 status rather than as a whole cohort,
## so this section works off explicit SAMPLE sets rather than cohort
## membership from here on - a stopifnot() partition check at the end
## verifies every kept sample landed in exactly one bucket.
##==========================
beta_to_mvalue <- function(beta, offset = 1e-6) {
  beta_clipped <- pmin(pmax(beta, offset), 1 - offset)
  log2(beta_clipped / (1 - beta_clipped))
}
mvalue_to_beta <- function(mvalue) {
  b <- 2^mvalue / (2^mvalue + 1)
  pmin(pmax(b, 0), 1)
}
 
drop_zero_var_probes <- function(mat, batch, label) {
  zero_var <- character(0)
  for (b in unique(batch)) {
    sub <- mat[, batch == b, drop = FALSE]
    v <- rowVars(sub); names(v) <- rownames(sub)
    zero_var <- union(zero_var, names(which(v == 0)))
  }
  if (length(zero_var) > 0) {
    cat(sprintf(
      "[%s] [%s] Dropping %d probes with zero variance within >=1 batch group\n",
      Sys.time(), label, length(zero_var)
    ))
    mat <- mat[!(rownames(mat) %in% zero_var), , drop = FALSE]
  }
  mat
}
 
run_combat <- function(mat, batch, label) {
  cat(sprintf("[%s] [%s] Running ComBat (%d probes x %d samples, batch levels: %s)...\n",
              Sys.time(), label, nrow(mat), ncol(mat), paste(unique(batch), collapse = ", ")))
  tryCatch({
    ComBat(dat = mat, batch = batch, mod = NULL,
           par.prior = TRUE, prior.plots = FALSE, BPPARAM = bpparam())
  }, error = function(e) {
    message(sprintf("[%s] ComBat with BPPARAM failed - retrying serially: %s",
                     label, conditionMessage(e)))
    ComBat(dat = mat, batch = batch, mod = NULL, par.prior = TRUE, prior.plots = FALSE)
  })
}
 
cohort_vector   <- samplesheet$Cohort[match(kept_samples, samplesheet$Sample_ID)]
platform_vector <- samplesheet$Platform[match(kept_samples, samplesheet$Sample_ID)]
stopifnot(identical(length(cohort_vector), ncol(final_matrix)))
stopifnot(!anyNA(cohort_vector), !anyNA(platform_vector))
 
known_cohorts <- c("Lucas", "MINT", "Capper", "Sturm")
unhandled_cohorts <- setdiff(unique(cohort_vector), known_cohorts)
if (length(unhandled_cohorts) > 0) {
  stop("Cohort(s) not routed by this section's logic: ",
       paste(unhandled_cohorts, collapse = ", "))
}
 
cat("\n[Cohort x Platform composition]\n")
print(table(cohort_vector, platform_vector))
 
cat(sprintf("[%s] Converting %d x %d beta matrix to M-values...\n",
            Sys.time(), nrow(final_matrix), ncol(final_matrix)))
mvalues <- beta_to_mvalue(final_matrix)
rm(final_matrix); gc(FALSE)
 
## ---- Identify the Sturm NF1-affected subset ------------------------------
## germline isn't carried into the merged samplesheet's columns, so this
## goes back to samplesheet_ref_pediatric directly (still in scope from
## Section 1). Intersected with kept_samples since the NF1 filter is
## applied to the FULL pediatric sheet, before Section 1's glioma-class
## filter and Pass 1's QC filter - some NF1 samples may not have survived
## either.
sturm_nf1_samples <- samplesheet_ref_pediatric %>%
  filter(grepl('Neurofibromatosis', germline)) %>%
  pull(Sample_Name)
sturm_nf1_samples <- intersect(sturm_nf1_samples, kept_samples[cohort_vector == "Sturm"])
 
cat(sprintf("\n[Sturm NF1] %d NF1-affected samples (of %d total kept Sturm samples)\n",
            length(sturm_nf1_samples), sum(cohort_vector == "Sturm")))
if (length(sturm_nf1_samples) < 3) {
  stop("Fewer than 3 Sturm NF1 samples survived filtering - ComBat cannot ",
       "reliably estimate a batch effect for this group in Step 2. Either ",
       "the NF1 filter matched too few samples, or most were dropped by an ",
       "earlier QC/class filter - check samplesheet_ref_pediatric$germline ",
       "and the glioma_pattern_ped/selected_classes filtering in Section 1.")
}
 
## ---- Step 1: within-MINT, batch = Platform ------------------------------
mint_samples  <- kept_samples[cohort_vector == "MINT"]
mint_platform <- platform_vector[cohort_vector == "MINT"]
cat(sprintf("\n[Step 1] Within-MINT platform correction: %d samples\n", length(mint_samples)))
print(table(mint_platform))
if (any(table(mint_platform) < 3)) {
  stop("A MINT platform group has <3 samples - ComBat cannot reliably ",
       "estimate a batch effect from that few samples.")
}
 
mvalues_mint <- mvalues[, mint_samples, drop = FALSE]
mvalues_mint <- drop_zero_var_probes(mvalues_mint, mint_platform, "Step1:MINT-platform")
mvalues_mint_corrected <- run_combat(mvalues_mint, mint_platform, "Step1:MINT-platform")
rm(mvalues_mint); gc(FALSE)
 
## ---- Step 2: Lucas + platform-corrected MINT + Sturm-NF1, batch=Cohort --
lucas_samples <- kept_samples[cohort_vector == "Lucas"]
cat(sprintf(
  "\n[Step 2] Lucas + MINT + Sturm-NF1 correction: %d Lucas + %d MINT + %d Sturm-NF1 = %d samples\n",
  length(lucas_samples), ncol(mvalues_mint_corrected), length(sturm_nf1_samples),
  length(lucas_samples) + ncol(mvalues_mint_corrected) + length(sturm_nf1_samples)
))
 
# All three pieces on the SAME probe subset that survived Step 1's
# zero-var drop (mvalues_mint_corrected's rows), so they align for cbind.
common_rows_step2 <- rownames(mvalues_mint_corrected)
mvalues_lucas     <- mvalues[common_rows_step2, lucas_samples, drop = FALSE]
mvalues_sturm_nf1 <- mvalues[common_rows_step2, sturm_nf1_samples, drop = FALSE]
 
mvalues_step2 <- cbind(mvalues_lucas, mvalues_mint_corrected, mvalues_sturm_nf1)
step2_batch <- c(
  rep("Lucas", length(lucas_samples)),
  rep("MINT", ncol(mvalues_mint_corrected)),
  rep("Sturm_NF1", length(sturm_nf1_samples))
)
if (any(table(step2_batch) < 3)) {
  stop("A Step 2 batch group has <3 samples - ComBat cannot reliably ",
       "estimate a batch effect from that few samples.")
}
 
mvalues_step2 <- drop_zero_var_probes(mvalues_step2, step2_batch, "Step2:Lucas-MINT-SturmNF1")
mvalues_step2_corrected <- run_combat(mvalues_step2, step2_batch, "Step2:Lucas-MINT-SturmNF1")
rm(mvalues_lucas, mvalues_sturm_nf1, mvalues_step2, mvalues_mint_corrected); gc(FALSE)
 
## ---- Uncorrected: Capper (all) + Sturm non-NF1 ---------------------------
final_probe_set <- rownames(mvalues_step2_corrected)  # probes surviving Step 1 + Step 2
capper_samples       <- kept_samples[cohort_vector == "Capper"]
sturm_all_samples    <- kept_samples[cohort_vector == "Sturm"]
sturm_nonnf1_samples <- setdiff(sturm_all_samples, sturm_nf1_samples)
uncorrected_samples  <- c(capper_samples, sturm_nonnf1_samples)
 
mvalues_uncorrected <- mvalues[final_probe_set, uncorrected_samples, drop = FALSE]
cat(sprintf(
  "\n[Uncorrected] %d samples (Capper: %d, Sturm non-NF1: %d) carried through UNCORRECTED\n",
  length(uncorrected_samples), length(capper_samples), length(sturm_nonnf1_samples)
))
 
## ---- Partition check: every kept sample accounted for exactly once ------
all_buckets <- c(mint_samples, lucas_samples, sturm_nf1_samples, uncorrected_samples)
stopifnot(
  length(all_buckets) == length(kept_samples),
  !any(duplicated(all_buckets)),
  setequal(all_buckets, kept_samples)
)
 
## ---- Combine into the final matrix, same column order as kept_samples ---
mvalues_combat <- cbind(mvalues_uncorrected, mvalues_step2_corrected)
mvalues_combat <- mvalues_combat[, kept_samples[kept_samples %in% colnames(mvalues_combat)], drop = FALSE]
rm(mvalues, mvalues_uncorrected, mvalues_step2_corrected); gc(FALSE)
 
cat(sprintf("[%s] Converting batch-corrected M-values back to beta values...\n", Sys.time()))
beta_matrix_combat <- mvalue_to_beta(mvalues_combat)
stopifnot(!anyNA(beta_matrix_combat))
 
## ---- Transparency log: what correction (if any) was applied per sample --
sample_order <- colnames(mvalues_combat)
correction_applied <- character(length(sample_order))
correction_applied[sample_order %in% capper_samples] <- "none (Capper - biologically distinct)"
correction_applied[sample_order %in% sturm_nonnf1_samples] <- "none (Sturm non-NF1 - biologically distinct)"
correction_applied[sample_order %in% mint_samples] <- "within-MINT (Platform) + Step2 (Lucas+MINT+SturmNF1)"
correction_applied[sample_order %in% lucas_samples] <- "Step2 (Lucas+MINT+SturmNF1)"
correction_applied[sample_order %in% sturm_nf1_samples] <- "Step2 (Lucas+MINT+SturmNF1)"
stopifnot(!any(correction_applied == ""))  # catches any sample this labeling missed
 
batch_correction_log <- data.frame(
  Sample_ID = sample_order,
  Cohort = cohort_vector[match(sample_order, kept_samples)],
  Platform = platform_vector[match(sample_order, kept_samples)],
  NF1_status = ifelse(sample_order %in% sturm_nf1_samples, "NF1", NA_character_),
  correction_applied = correction_applied,
  row.names = NULL
)
fwrite(batch_correction_log, file.path(LOG_DIR, "batch_correction_log.csv"))
 
##==========================
## 12. SAVE BATCH-CORRECTED OUTPUTS
##==========================
saveRDS(mvalues_combat, file.path(OUT_DIR, "mvalue_matrix_ComBat_nested.rds"))
saveRDS(beta_matrix_combat, file.path(OUT_DIR, "beta_matrix_ComBat_nested.rds"))
if (WRITE_CSV_EXPORTS) {
  fwrite(
    data.table(Probe_ID = rownames(beta_matrix_combat), beta_matrix_combat),
    file.path(OUT_DIR, "beta_matrix_ComBat_nested.csv.gz")
  )
}
 
cat(sprintf(
  "[%s] Batch-corrected beta matrix saved: %d probes x %d samples\n",
  Sys.time(), nrow(beta_matrix_combat), ncol(beta_matrix_combat)
))
cat("  - logs/batch_correction_log.csv (per-sample: what correction was applied)\n")
 

##==========================
## 12b. WRITE ANNDATA (.h5ad)
##==========================
## Uncorrected betas in X, batch-corrected matrices as layers.
##
## The ComBat output covers FEWER probes than the uncorrected matrix (Step 1 and
## Step 2 each drop probes with zero variance within a batch group), so the
## layers are padded back to the full `complete_probes` space with NA rather than
## shrinking X to the intersection. var$in_combat marks which probes are real.
## Every kept sample is present in both, guaranteed by the Section 11 partition
## check, so the columns need no padding.
final_matrix <- readRDS(file.path(OUT_DIR, "beta_matrix_complete_HM450space.rds"))
final_matrix <- final_matrix[, kept_samples, drop = FALSE]

pad_to_full <- function(m) {
  out <- matrix(NA_real_, nrow = length(complete_probes), ncol = length(kept_samples),
                dimnames = list(complete_probes, kept_samples))
  out[rownames(m), colnames(m)] <- m
  out
}
## Only the M-value layer is stored; beta is recoverable via mvalue_to_beta()
## and a second full-size matrix is expensive at this scale.
mvalue_combat_full <- pad_to_full(mvalues_combat)

cat(sprintf("[%s] ComBat covers %d/%d probes; the remainder are NA in the layers\n",
            Sys.time(), nrow(beta_matrix_combat), length(complete_probes)))

## ---- obs: one row per kept sample, in matrix column order ----------------
## rownames are set explicitly: without them pandas assigns a positional
## RangeIndex and the sample identifiers survive only as a column.
obs <- samplesheet %>%
  filter(Sample_ID %in% kept_samples) %>%
  as.data.frame()
obs <- obs[match(kept_samples, obs$Sample_ID), , drop = FALSE]
## Columns are selected explicitly rather than joining whole frames: sample_qc
## also carries Platform, which would collide with obs$Platform and silently
## become Platform.x / Platform.y.
##
## frac_missing_* are NA-fractions from getBetas(mask = FALSE), i.e. genuinely
## unmeasurable signal (plus, for _hm450, structural mLiftOver gaps) - NOT
## pOOBAH detection failure. frac_failed_detection is the detection-based rate
## from the Section 5b scan, on native probes.
obs <- obs %>%
  left_join(batch_correction_log %>% select(Sample_ID, correction_applied),
            by = "Sample_ID")  %>%
    left_join(sample_detection_qc %>% group_by(Sample_ID) %>% filter(row_number() == 1) %>% ungroup(), by = "Sample_ID")


##==========================
## 12a. SAMPLE ANNOTATION: classes, families, pathology, location, age clock
##==========================
## Lifted verbatim out of Create_embedding_df() in the original Section 14, so
## that every annotation lives in obs rather than being recomputed inside each
## embedding call. The case_when() chains, the join_key construction and the
## Classes lookup are unchanged; only their location has moved.
##
## Horvath age is estimated on the UNCORRECTED betas, matching the original
## (which read beta_matrix_complete_HM450space.rds for this).
beta_for_clock <- final_matrix
age_estimates <- methylclock::DNAmAge(beta_for_clock, min.perc = 0.6, clocks = 'Horvath')
rm(beta_for_clock); gc(FALSE)

Classes <- read.csv(classes_csv) %>%
    rename(Methylation_class_MINT = Methylation_Class) %>% select(-Cohort) %>%
    filter(Sample_ID != '') %>% select(Sample_ID,Loc,PA,Methylation_class_MINT) %>%
    rbind(samplesheet_ref_pediatric %>% mutate(PA = who_type,Sample_ID = Sample_Name, Methylation_class_MINT = meth_class_rf, Loc = tumor_location)  %>%
          select(Sample_ID,Loc,PA,Methylation_class_MINT))  %>%
    unique()

obs <- obs %>%
        mutate(
            m_num = str_extract(Sample_ID, "(?<=_M)\\d+"),
            r_num = str_extract(Sample_ID, "(?<=_R)\\d+"),
            join_key =ifelse(Cohort == 'MINT', paste0("M_", m_num, ifelse(r_num == "1", "", paste0("R", r_num))),Sample_ID)) %>%
    left_join(Classes, by = c("join_key" = "Sample_ID")) %>%
    left_join(age_estimates, by = c('Sample_ID' = 'id')) %>% 
    mutate(
        Methylation_Class = ifelse(!is.na(Methylation_class_MINT),Methylation_class_MINT,Methylation_Class),
        family = case_when(
            is.na(Methylation_Class) ~ 'Unclass.',
            Methylation_Class == 'PA, NF1-associated' ~ 'Unclass.',
            grepl('^CONTR_|^Control,', Methylation_Class, ignore.case = TRUE) ~ 'Control',
            grepl('^GBM|^High-grade glioma|^Posterior fossa glioblastoma', 
                  Methylation_Class, ignore.case = TRUE) ~ 'GBM',
            grepl('^DMG|^Diffuse midline glioma', Methylation_Class, ignore.case = TRUE) ~ 'DMG',
            grepl('^O_IDH$|^A_IDH|^Oligodendroglioma, IDH-mutated|^Astrocytoma, IDH-mutated|^High-grade astrocytoma, IDH-mutated', 
                  Methylation_Class, ignore.case = TRUE) ~ 'IDHmt',
            Methylation_Class == 'HGAP' | 
            grepl('ANA_PA|^Anaplastic pilocytic astrocytoma', Methylation_Class, ignore.case = TRUE) ~ 'HGAP',
            grepl('^LGG_PA', Methylation_Class) ~ 'LGG PA',
            grepl('^PA|^Pilocytic Astrocytoma', Methylation_Class, ignore.case = TRUE) ~ 'PA',
            Methylation_Class == 'LGG_RGNT' | 
            grepl('^Rosette-forming glioneuronal', Methylation_Class, ignore.case = TRUE) ~ 'RGNT',
            Methylation_Class == 'PXA' | 
            grepl('^Pleomorphic xanthoastrocytoma', Methylation_Class, ignore.case = TRUE) ~ 'PXA',
            grepl('^LGG|^Dysembryoplastic neuroepithelial|^Subependymal giant cell astrocytoma|^Desmoplastic infantile ganglioglioma|^Ganglioglioma$|^Low-grade glioma', 
                  Methylation_Class, ignore.case = TRUE) ~ 'LGG',
            grepl('^IHG$|^Infantile hemispheric glioma', Methylation_Class, ignore.case = TRUE) ~ 'IHG',
            grepl('^PGG', Methylation_Class) ~ 'PGG',
            Methylation_Class == 'Non-classifiable' ~ 'Unclass.',
            TRUE ~ 'Unclass.'   # catches "Non-classifiable" and anything unmapped
        ),
        
        Subtype_Pathology = case_when(
            is.na(PA) ~ 'Unclass.',
            PA %in% c('Descriptive diagnosis', 'Non-classifiable') ~ 'Unclass.',
            PA == 'Non-neoplastic tissue'  ~ 'Control',
            
            grepl('anaplastic pilocytic astrocytoma', PA, ignore.case = TRUE) ~ 'HGAP',
            grepl('high-grade astrocytoma with piloid features', PA, ignore.case = TRUE) ~ 'HGAP',
            
            grepl('pilocytic astrocytoma|pilomyxoid astrocytoma', PA, ignore.case = TRUE) ~ 'PA',
            
            grepl('rosette-forming glioneuronal', PA, ignore.case = TRUE) ~ 'RGNT',
            grepl('multinodular and vacuolating neuronal tumor', PA, ignore.case = TRUE) ~ 'MVNT',
            grepl('xanthoastrocytoma', PA, ignore.case = TRUE) ~ 'PXA',
            grepl('subependymal giant cell astrocytoma', PA, ignore.case = TRUE) ~ 'SEGA',
            grepl('^ganglioglioma$|^anaplastic ganglioglioma$', PA, ignore.case = TRUE) ~ 'GGL',
            PA == 'Diffuse glial tumor' ~ 'Diffuse glioma, NOS',
            
            grepl('diffuse midline glioma', PA, ignore.case = TRUE) ~ 'DMG',
            
            grepl('glioblastoma.*idh-wildtype|glioblastoma \\(gbm\\)|glioblastoma, h3 g34|gliosarcoma', 
                  PA, ignore.case = TRUE) ~ 'GBM',
            
            grepl('idh-mutant|idh-mutated|1p/19q-codeleted|^oligoastrocytoma$|^anaplastic oligoastrocytoma$', 
                  PA, ignore.case = TRUE) ~ 'IDHmt',
            
            grepl('^high-grade astrocytoma$|^high-grade glioma, nos$|^anaplastic astrocytoma$|
             ^anaplastic astrocytoma\\. idh-wildtype$|^anaplastic oligodendroglioma$', 
             PA, ignore.case = TRUE) ~ 'HGG, NOS',
            
            grepl('dysembryoplastic neuroepithelial|desmoplastic infantile|low-grade astrocytoma|
             paediatric diffuse astrocytoma|angiocentric|polymorphous low-grade neuroepithelial|
             ^diffuse astrocytoma$', PA, ignore.case = TRUE) ~ 'LGG',
            
            grepl('papillary glioneuronal tumor', PA, ignore.case = TRUE) ~ 'PGNT',
            
            grepl('medulloblastoma', PA, ignore.case = TRUE) ~ 'Medulloblastoma',
            grepl('^ependymoma', PA, ignore.case = TRUE) ~ 'Ependymoma',
            grepl('meningioma', PA, ignore.case = TRUE) ~ 'Meningioma',
            grepl('craniopharyngioma', PA, ignore.case = TRUE) ~ 'Craniopharyngioma',
            grepl('germinoma|germ cell tumor|teratoma', PA, ignore.case = TRUE) ~ 'Germ cell tumor',
            grepl('pineoblastoma|pineocytoma|pineal parenchymal tumor', PA, ignore.case = TRUE) ~ 'Pineal tumor',
            grepl('embryonal tumor with multilayered rosettes|cns embryonal tumour|
             primitive neuroectodermal tumor|cns neuroblastoma|ganglioneuroblastoma', 
             PA, ignore.case = TRUE) ~ 'Embryonal tumor',
            
            TRUE ~ 'Other'
        ),
        Loc = case_when(
      is.na(Loc) ~ 'Unknown',
      
      # --- Spinal: check before supratentorial/posterior fossa since some rows combine multiple compartments ---
      grepl('^spinal|spinal cord', Loc, ignore.case = TRUE) ~ 'Spinal',
      
      # --- Posterior fossa / infratentorial (brainstem, cerebellum, pons, medulla, IV ventricle) ---
      grepl('posterior fossa|infratentorial|cerebellum|brainstem|^pons$|medulla oblongata|
             fourth ventricle|ventricle iv|cerebellar', Loc, ignore.case = TRUE) ~ 'Posterior fossa',
      
      # --- Supratentorial (covers hemispheric, diencephalic, sellar, intraventricular II/III, and the short labels) ---
      grepl('supratentorial|hemispheric|diencephalic|thalam|hypothalam|frontal|temporal|parietal|
             occipital|sellar|intraventricular|ventricle ii|ventricle iii|basal ganglia|
             basal nuclei|corpus callosum|pineal region|orbital frontal|somatosensory|motor cortex', 
            Loc, ignore.case = TRUE) ~ 'Supratentorial',
      
      # --- Peripheral / cranial nerve, meningeal, other non-parenchymal sites ---
      grepl('peripheral nerve|cranial nerve|meningeal|mesenchyme|cranium|dermal', 
            Loc, ignore.case = TRUE) ~ 'Extra-axial/Peripheral',
      
      TRUE ~ 'Other/Unclassified'
    )


        
    )


## ---- NF1 status ----------------------------------------------------------
## Cohort-level assignment, except for Sturm which splits per sample on
## constitutional NF1 status (sturm_nf1_samples, derived in Section 11 from
## samplesheet_ref_pediatric$germline). NF1_status_source records how each
## label was arrived at, since only the Sturm calls rest on a per-sample
## germline annotation - the rest are assumptions about the cohort as a whole.
obs <- obs %>%
  mutate(
    NF1_status = case_when(
      Cohort %in% c("Lucas", "MINT")                        ~ "NF1",
      Cohort == "Capper"                                    ~ "NF1wt",
      Cohort == "Sturm" & Sample_ID %in% sturm_nf1_samples  ~ "NF1",
      Cohort == "Sturm"                                     ~ "NF1wt",
      TRUE                                                  ~ NA_character_
    ),
    NF1_status_source = case_when(
      Cohort == "Sturm" ~ "germline annotation (per sample)",
      TRUE              ~ "cohort assumption"
    )
  )
stopifnot(!any(is.na(obs$NF1_status)))

cat("\n[NF1 status by cohort]\n")
print(table(obs$Cohort, obs$NF1_status))

## Capper is labelled NF1wt wholesale, but the reference set does contain an
## NF1-associated pilocytic astrocytoma methylation class - flag any such
## samples rather than let the cohort-level assumption bury them.
capper_nf1_class <- obs$Sample_ID[
  obs$Cohort == "Capper" & grepl("NF1", obs$Methylation_Class, ignore.case = TRUE)
]
if (length(capper_nf1_class) > 0) {
  warning(sprintf(
    paste0("%d Capper sample(s) carry an NF1-associated methylation class but ",
           "are labelled NF1wt by the cohort-level rule: %s"),
    length(capper_nf1_class),
    paste(head(capper_nf1_class, 5), collapse = ", ")
  ))
}

## Everything Create_embedding_df() used to derive is now a column in obs.
cat("\n[Annotation summary: family]\n")
print(table(obs$family, useNA = "ifany"))
cat("\n[Annotation summary: Subtype_Pathology]\n")
print(table(obs$Subtype_Pathology, useNA = "ifany"))
cat("\n[Annotation summary: Loc]\n")
print(table(obs$Loc, useNA = "ifany"))

rownames(obs) <- obs$Sample_ID
stopifnot(identical(rownames(obs), colnames(final_matrix)))

## ---- var: one row per probe ---------------------------------------------
var <- data.frame(
  Probe_ID  = complete_probes,
  in_combat = complete_probes %in% rownames(beta_matrix_combat),
  row.names = complete_probes,
  stringsAsFactors = FALSE
)

## ---- Assemble and write --------------------------------------------------
## AnnData is samples x probes, so each matrix is transposed on the way in.
adata <- AnnData(
  X      = t(final_matrix),
  obs    = obs,
  var    = var,
  layers = list(
    mvalue_combat = t(mvalue_combat_full)
  ),
  uns    = list(
    kept_samples        = kept_samples,
    complete_probes     = complete_probes,
    combat_probes       = rownames(beta_matrix_combat),
    sturm_nf1_samples   = sturm_nf1_samples,
    mint_samples        = mint_samples,
    lucas_samples       = lucas_samples,
    uncorrected_samples = uncorrected_samples,
    retention_report    = retention_df,
    target_platform     = TARGET_PLATFORM,
    max_missing_fraction = MAX_MISSING_FRACTION,
    drop_sex_chrom      = DROP_SEX_CHROM
  )
)

dir.create(dirname(adata_out), recursive = TRUE, showWarnings = FALSE)
write_h5ad(adata, adata_out)

cat(sprintf("[%s] Wrote %s: %d samples x %d probes (X = uncorrected beta; layers = mvalue_combat)\n",
            Sys.time(), adata_out, nrow(obs), nrow(var)))

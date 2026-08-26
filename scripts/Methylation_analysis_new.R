#!/usr/bin/env Rscript
## =============================================================================
## MULTI-PLATFORM METHYLATION HARMONIZATION PIPELINE
## sesame + mLiftOver | 450K / EPICv1 / EPICv2 -> common HM450 probe space
## =============================================================================

##==========================
## 0. PACKAGES, PARALLEL BACKEND, PATHS
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
library(data.table)
library(methylclock)
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 16))
register(MulticoreParam(workers = ncores, stop.on.error = FALSE), default = TRUE)

if (exists("snakemake")) {
  idat_dir_ref       <- snakemake@params[['idat_dir_Capper']]
  idat_dir_lucas     <- snakemake@params[['idat_dir_Lucas']]
  samplesheet_ref    <- snakemake@params[['samplesheet_Capper']]
  samplesheet_lucas  <- snakemake@params[['samplesheet_Lucas']]
  Zhou_input         <- snakemake@params[['Zhou_probes']]
  CrossReactive_input<- snakemake@params[['CrossReactive_probes']]
  Problematic_input  <- snakemake@params[['Problematic_probes']]
  samplesheet_mint    <- snakemake@input[[1]]
  OUT_DIR             <- snakemake@params[['out_dir']]
} else {
    idat_dir_ref     <- '~/mnt/BIGR_home/MINT/data/idat/Capper'
    
    
  idat_dir_lucas   <- '~/mnt/BIGR_home/MINT/data/idat/Lucas'
  idat_dir_mint    <- '~/mnt/BIGR_home/MINT/data/idat/'
  samplesheet_mint <- '~/mnt/BIGR_home/MINT/output/Methylation/samplesheets/Samplesheet_Methylation.csv'
    samplesheet_ref  <- '~/mnt/BIGR_home/MINT/data/samplesheet_capper.csv'
    samplesheet_ref_pediatric  <- '~/mnt/BIGR_home/MINT/data/Samplesheet_Sturm.csv'

  samplesheet_lucas<- '~/mnt/BIGR_home/MINT/data/SampleSheetLucas.csv'

  Zhou_input          <- '~/mnt/BIGR_home/Resources/EPIC/manifest/AppendixD_Zhou_et_al_MASKgeneral_list.txt'
  CrossReactive_input <- '~/mnt/BIGR_home/Resources/EPIC/manifest/AppendixE_CrossReactiveProbes_EPICv1.txt'
  Problematic_input   <- '~/mnt/BIGR_home/Resources/EPIC/manifest/AppendixF_ProblematicProbes_EPICv1-b5.txt'

  OUT_DIR <- '~/mnt/BIGR_home/MINT/output/Methylation/harmonized'
}

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
  rowSums(failed_hm450, na.rm = TRUE)
}

n_failed_detection_per_probe <- integer(length(hm450_probe_ids))
names(n_failed_detection_per_probe) <- hm450_probe_ids

for (plat in c("HM450", "EPIC", "EPICv2")) {
  plat_sheet <- samplesheet %>% filter(Platform == plat, Sample_ID %in% kept_samples)
  if (nrow(plat_sheet) == 0) next
  n_chunks <- ceiling(nrow(plat_sheet) / CHUNK_SIZE)
  cat(sprintf("[%s] Detection-p-value scan %s: %d samples in %d chunk(s)\n",
              Sys.time(), plat, nrow(plat_sheet), n_chunks))
  for (i in seq_len(n_chunks)) {
    idx <- ((i - 1) * CHUNK_SIZE + 1):min(i * CHUNK_SIZE, nrow(plat_sheet))
    chunk_failed_counts <- compute_detection_failcount_chunk(
      prefixes   = plat_sheet$Basename_full[idx],
      sample_ids = plat_sheet$Sample_ID[idx],
      platform   = plat
    )
    n_failed_detection_per_probe <- n_failed_detection_per_probe + chunk_failed_counts
    gc(FALSE)
  }
}
frac_failed_detection_per_probe <- n_failed_detection_per_probe / length(kept_samples)
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
## 13. DIAGNOSTIC: PCA before/after, colored by Cohort
##==========================
diagnostic_pca <- function(mat, cohort, label, top_n = 5000) {
  vars <- rowVars(mat)
  names(vars) <- rownames(mat)
  top_probes <- names(sort(vars, decreasing = TRUE))[seq_len(min(top_n, length(vars)))]
  pca <- prcomp(t(mat[top_probes, ]), scale. = FALSE)
  df <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], Cohort = cohort)
  ggplot(df, aes(PC1, PC2, color = Cohort)) +
    geom_point(alpha = 0.6, size = 1.5) +
    labs(title = label) +
    theme_minimal()
}
 
beta_precorrect <- readRDS(file.path(OUT_DIR, "beta_matrix_complete_HM450space.rds"))
common_probes_for_plot <- intersect(rownames(beta_precorrect), rownames(mvalues_combat))
p_before <- diagnostic_pca(
  beta_to_mvalue(beta_precorrect[common_probes_for_plot, ]),
  cohort_vector, "M-values BEFORE correction"
)
p_after <- diagnostic_pca(
  mvalues_combat[common_probes_for_plot, ],
  cohort_vector, "M-values AFTER correction (Capper + Sturm non-NF1 uncorrected)"
)
rm(beta_precorrect); gc(FALSE)
 
ggsave(file.path(LOG_DIR, "PCA_before_ComBat.png"), p_before, width = 6, height = 5)
ggsave(file.path(LOG_DIR, "PCA_after_ComBat.png"),  p_after,  width = 6, height = 5)
cat("  - logs/PCA_before_ComBat.png / PCA_after_ComBat.png (visual sanity check)\n")

find_elbow <- function(y) {
  # Create an index vector for the x-axis
  x <- 1:length(y)
  
  # Coordinates of the first and last points
  p1 <- c(x[1], y[1])
  p2 <- c(x[length(x)], y[length(y)])
  
  # Vector defining the line between endpoints
  line_vec <- p2 - p1
  
  # Calculate squared orthogonal distance for each point
  distances <- sapply(x, function(i) {
    point <- c(x[i], y[i])
    p1_to_point <- point - p1
    
    # Vector projection math
    projection <- sum(p1_to_point * line_vec) / sum(line_vec^2)
    orthogonal_vec <- p1_to_point - projection * line_vec
    
    return(sum(orthogonal_vec^2)) # Squared distance is enough to find the max
  })
  
  # Return the index of the maximum distance
  return(which.max(distances))
}




##==========================
## 14. DIMENSIONALITY REDUCTION FOR VISUALIZATION: top-variable probes ->
##     PCA (100 PCs) -> t-SNE / UMAP
##==========================
## Uses the FINAL, batch-corrected beta matrix (beta_matrix_combat) - the
## endpoint of the pipeline. If you'd rather use the pre-ComBat betas
## instead, swap in: readRDS(file.path(OUT_DIR, "beta_matrix_complete_HM450space.rds"))
##
## Requires the `uwot` package for UMAP (not loaded above) - install with
## install.packages("uwot") if missing. PCA uses base-R prcomp() - computes
## a full SVD (up to min(n_samples-1, n_probes) components) then keeps only
## the first N_PCS. Fine at this scale (~2864 samples), just note this is
## more compute than strictly necessary for only wanting 100 components -
## an irlba-based truncated SVD would be faster if this step becomes slow.
if (!requireNamespace("uwot", quietly = TRUE)) {
  stop("Package 'uwot' is required for UMAP - install.packages('uwot')")
}
library(uwot)

Run_PCA <- function(input_data){
    beta_matrix_combat <- mvalue_to_beta(input_data)
    probe_vars <- rowVars(beta_matrix_combat)
    N_TOP_VARIABLE_PROBES <- find_elbow(probe_vars)
    print(paste0('Number of variable probes: ',N_TOP_VARIABLE_PROBES))
    names(probe_vars) <- rownames(beta_matrix_combat)
    top_var_probes <- names(sort(probe_vars, decreasing = TRUE))[
        seq_len(min(N_TOP_VARIABLE_PROBES, length(probe_vars)))
    ]
    mat_top <- beta_matrix_combat[top_var_probes, , drop = FALSE]
    pca_full <- prcomp(t(mat_top), center = TRUE, scale. = FALSE)
    return(pca_full)
}

beta_for_clock <- readRDS(file.path(OUT_DIR, "beta_matrix_complete_HM450space.rds"))
age_estimates <- methylclock::DNAmAge(beta_for_clock, min.perc = 0.6,clocks = 'Horvath')
rm(beta_for_clock)

Create_embedding_df<- function(input_data, perplexity = 30, N_PCS = 100){
    pca_full <- Run_PCA(input_data)
    pca_scores <- pca_full$x[, seq_len(min(N_PCS, ncol(pca_full$x))), drop = FALSE]
    set.seed(42)
    tsne_result <- Rtsne(
        pca_scores, dims = 2, pca = FALSE, check_duplicates = FALSE,
        perplexity = min(perplexity, floor((nrow(pca_scores) - 1) / 3)),
        verbose = TRUE, num_threads = ncores
    )
    tsne_coords <- tsne_result$Y
    rownames(tsne_coords) <- rownames(pca_scores)
    colnames(tsne_coords) <- c("tSNE1", "tSNE2")

    Classes <- read.csv('/home/jurriaan/mnt/BIGR_home/MINT/data/SamplesheetCapper_with_classes.csv') %>%
        rename(Methylation_class_MINT = Methylation_Class) %>% select(-Cohort) %>%
        filter(Sample_ID != '') %>% select(Sample_ID,Loc,PA,Methylation_class_MINT) %>%
        rbind(samplesheet_ref_pediatric %>% mutate(PA = who_type,Sample_ID = Sample_Name, Methylation_class_MINT = meth_class_rf, Loc = tumor_location)  %>%
              select(Sample_ID,Loc,PA,Methylation_class_MINT))  %>%
        unique() 

    ## ---- Combine + save --------------------------------------------------
    embedding_df <- data.frame(
        Sample_ID = rownames(pca_scores),
        PC1 = pca_scores[,1],
        PC2 = pca_scores[,2],
        PC3 = pca_scores[,3],
        PC4 = pca_scores[,4],

        Cohort = samplesheet$Cohort[match(rownames(pca_scores), samplesheet$Sample_ID)],
        Platform = samplesheet$Platform[match(rownames(pca_scores), samplesheet$Sample_ID)],
        Methylation_Class = samplesheet$Methylation_Class[match(rownames(pca_scores), samplesheet$Sample_ID)],
        Mean_methylation = colMeans(input_data),
        tSNE1 = tsne_coords[, 1], tSNE2 = tsne_coords[, 2],
        #UMAP1 = umap_coords[, 1], UMAP2 = umap_coords[, 2],
        
        row.names = NULL
    ) %>%
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

    
    return(embedding_df)
}


## ---- Plots -------------------------------------------------------------
color_pal <- c("#4E79A7","#A0CBE8","#F28E2B","#FFBE7D","#59A14F","#8CD17D","#B6992D","#F1CE63","#499894","#86BCB6","#E15759","#FF9D9A","#79706E","#BAB0AC","#D37295","#FABFD2","#B07AA1","#D4A6C8","#9D7660","#D7B5A6")
names(color_pal) <- c("Blue","Light Blue" ,"Orange","Light Orange","Green","Light Green","Yellow-Green","Yellow","Teal","Light Teal", "Red","Pink","Dark Gray","Light Gray","Pink","Light Pink","Purple","Light Purple","Brown","Light Orange")
cohort_pal <- list('Capper' = color_pal[['Light Gray']],'Lucas' = color_pal[['Red']],'MINT' = color_pal[['Green']], 'Sturm' = color_pal[['Orange']], 'Sturm NF1' = color_pal[['Orange']], 'Sturm NF1wt' = color_pal[['Yellow']])
loc_pal <- c("#4E79A7","#F28E2B","#76B7B2","#59A14F","#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC")
subtype_pal <- list('GBM' = color_pal[['Yellow-Green']], DMG = color_pal[['Brown']], 'PA' = color_pal[['Blue']],
                    'RGNT' = color_pal[['Teal']], 'LGG' = color_pal[['Light Teal']], 'IDHmt' = color_pal[['Light Orange']], 'PXA' = color_pal[['Purple']], 'HGAP' = color_pal[['Yellow']],
                    'PGG' = color_pal[['Pink']] , 'IHG' = color_pal[['Red']],'Control' = color_pal[['Dark Gray']], 'Unclass.' = 'white', 'PA, NF1' = color_pal[['Green']])

save.image(file = "AllData.RData")
load("dump/AllData.RData")




beta_values_combat <- mvalue_to_beta(mvalues_combat)
embeddings_all <- Create_embedding_df(mvalues_combat)

write.csv(embeddings_all %>%
            select(Sample_ID,Cohort,Platform,tSNE1,tSNE2,Loc,family,Subtype_Pathology),
            file = 'Embeddings_and_subtypes.csv',row.names = F, quote=F)


Sample_selection <- c(
    samplesheet_ref_pediatric %>% filter(grepl('Neurofibromatosis', germline)) %>% pull(Sample_Name),
    embeddings_all %>% filter(Cohort %in% c('Capper','Lucas','MINT')) %>% pull(Sample_ID)
)

embeddings_all_selected <- Create_embedding_df(mvalues_combat[,Sample_selection]) %>%
    mutate(Cohort = ifelse(Cohort == 'Sturm','Sturm NF1', Cohort)) %>%
    arrange(Cohort)

write.csv(embeddings_all_selected %>%
            select(Sample_ID,Cohort,Platform,tSNE1,tSNE2,Loc,family,Subtype_Pathology),
            file = 'Embeddings_and_subtypes_selected.csv',row.names = F, quote=F)

tsne_family <- embeddings_all %>%
  ggplot(aes(tSNE1, tSNE2, fill = family, shape = Cohort, color = Cohort)) +
  geom_point(alpha = 0.9, size = 3,stroke = 1) +
    scale_fill_manual(values = subtype_pal) +
    scale_color_manual(values = cohort_pal) +
  scale_shape_manual(values = c(22,24,21,23)) +
  theme_classic(base_size = 13) +
  labs(
       fill = 'DKFZ class') +
  guides(
    fill = guide_legend(override.aes = list(shape = 21)),   # ensures fill legend uses a fillable shape
    shape = guide_legend(override.aes = list(fill = "white"))  # gives shape legend keys a visible fill
  )  
    
pdf('tSNE_AllGliomas_w_batch_correction.pdf', height = 6 , width = 8)
tsne_family
dev.off()



tsne_family_zoom <- embeddings_all_selected %>%
    filter(tSNE1 > -20, tSNE1 < 0,tSNE2 > -20, tSNE2 < 0 ) %>%
  ggplot(aes(tSNE1, tSNE2, fill = family, shape = Cohort, color = Cohort)) +
  geom_point(alpha = 0.9, size = 4,stroke = 1) +
    scale_fill_manual(values = subtype_pal) +
    scale_color_manual(values = cohort_pal) +
  scale_shape_manual(values = c(22,24,21,23)) +
  theme_classic(base_size = 13) +
  labs(
       fill = 'DKFZ class') +
  guides(
    fill = guide_legend(override.aes = list(shape = 21)),   # ensures fill legend uses a fillable shape
    shape = guide_legend(override.aes = list(fill = "white"))  # gives shape legend keys a visible fill
  ) +
    theme(legend.position = 'none')
    
pdf('tSNE_AllGliomas_w_batch_correction_zoom_tsne_-20_0.pdf', height = 5 , width = 5)
tsne_family_zoom
dev.off()



tsne_Subtype_Pathology <- embeddings_all_selected %>%
  ggplot(aes(tSNE1, tSNE2, fill = Subtype_Pathology, shape = Cohort, color = Cohort)) +
  geom_point(alpha = 0.9, size = 2,stroke = 1) +
    scale_fill_manual(values = subtype_pal) +
    scale_color_manual(values = cohort_pal) +
  scale_shape_manual(values = c(22,24,21,23)) +
  theme_classic(base_size = 13) +
  labs(
       fill = 'DKFZ class') +
  guides(
    fill = guide_legend(override.aes = list(shape = 21)),   # ensures fill legend uses a fillable shape
    shape = guide_legend(override.aes = list(fill = "white"))  # gives shape legend keys a visible fill
  ) 

pdf('tSNE_AllGliomas_w_batch_correction_Histological_subtype.pdf', height = 6 , width = 8)
tsne_Subtype_Pathology
dev.off()

PA_samples <- embeddings_all_selected %>%
    filter(Cohort %in% c('MINT','Lucas','Sturm NF1'),
           family %in% c('PA','PA, NF1','Control') |
           (family == 'Unclass.' & Subtype_Pathology == 'PA') |
           (family == 'Control' & Subtype_Pathology == 'PA')) %>%
    pull(Sample_ID) %>% unique()

embeddings_PA <- Create_embedding_df(mvalues_combat[,PA_samples] , perplexity = 10) %>%
    mutate(Cohort = ifelse(Cohort == 'Sturm','Sturm NF1', Cohort)) %>%
    arrange(Cohort)
pca_full_PA<- Run_PCA(mvalues_combat[,PA_samples])
var_explained <- (pca_full_PA$sdev^2) / sum(pca_full_PA$sdev^2) * 100

pdf('PCA_PA_gliomas_w_batch_correction_family.pdf', height = 3, width = 6)
PCA_PA <- embeddings_PA %>%
    ggplot(aes(PC1, PC2, fill = family, shape = Cohort, color = Cohort)) +
    geom_point( size = 3,stroke = 1) +
    scale_fill_manual(values = subtype_pal) +
    scale_color_manual(values = cohort_pal) +
    scale_shape_manual(values = c(22,24,21,23)) +
    theme_classic(base_size = 13) +
    labs(
        fill = 'DKFZ class',
        x = paste0('PC1 (',round(var_explained[1]),'%)' ),
        y = paste0('PC2 (',round(var_explained[2]),'%)' ))+
    guides(
        fill = guide_legend(override.aes = list(shape = 21))) +
    geom_point(data = embeddings_PA %>% filter(Methylation_Class == 'PA, NF1-associated'),aes(PC1,PC2), shape = 8, color = 'black', inherit.aes = F, size = 1)
PCA_PA
dev.off()


pdf('PCA_PA_gliomas_w_batch_correction_Location.pdf', height = 3, width = 6)
PCA_PA_Loc <- embeddings_PA %>%
    ggplot(aes(PC1, PC2, fill = Loc, shape = Cohort, color = Cohort)) +
    geom_point( size = 3,stroke = 1) +
    scale_fill_manual(values = loc_pal) +
    scale_color_manual(values = cohort_pal) +
    scale_shape_manual(values = c(22,24,21,23)) +
    theme_classic(base_size = 13) +
    labs(
        fill = 'DKFZ class',
        x = paste0('PC1 (',round(var_explained[1]),'%)' ),
        y = paste0('PC2 (',round(var_explained[2]),'%)' ))+
    guides(
        fill = guide_legend(override.aes = list(shape = 21))) +
    geom_point(data = embeddings_PA %>% filter(Methylation_Class == 'PA, NF1-associated'),aes(PC1,PC2), shape = 8, color = 'black', inherit.aes = F, size = 1)
PCA_PA_Loc
dev.off()



write.csv(embeddings_PA %>%
            select(Sample_ID,Cohort,Platform,PC1,PC2,Loc,family,Subtype_Pathology),
            file = 'PA_PCA_and_subtypes.csv',row.names = F, quote=F)


## =============================================================================
## In-silico PURIFICATION (subtraction), plotted the SAME way as the
## additive dilution series - combined arrow plot, faceted by fraction,
## colored by reference cell type, using the SAME ratios (0.25/0.5/0.75/1).
##
## Depends on make_purification_series() and ref_signatures already being
## defined (see PA_purification_analysis.R). This script only changes the
## ratio grid and the plotting style to match your dilution series - the
## underlying subtraction formula (purified = (observed - f*ref)/(1-f)) is
## unchanged.
## =============================================================================
make_purification_series <- function(tumor_sample_id, ref_signatures, cell_type,
                                      beta_matrix, pca_object,
                                      fractions = seq(0, 0.98, by = 0.02)) {
  probe_set <- rownames(pca_object$rotation)
  stopifnot(tumor_sample_id %in% colnames(beta_matrix))
  stopifnot(cell_type %in% colnames(ref_signatures))
 
  common_probes <- Reduce(intersect, list(
    probe_set, rownames(beta_matrix), rownames(ref_signatures)
  ))
  if (length(common_probes) < length(probe_set) * 0.5) {
    warning(sprintf(
      "[%s/%s] Only %d/%d PCA probes covered by the reference - purification ",
      "may be unreliable for this cell type.",
      tumor_sample_id, cell_type, length(common_probes), length(probe_set)
    ))
  }
 
  beta_tumor <- beta_matrix[common_probes, tumor_sample_id]
  ref_vec    <- ref_signatures[common_probes, cell_type]
 
  # Standard linear demixing: purified = (observed - f*ref) / (1-f)
  purified_betas <- sapply(fractions, function(f) {
    p <- (beta_tumor - f * ref_vec) / (1 - f)
    pmin(pmax(p, 0), 1)  # clip to valid beta range - demixing can overshoot
                         # [0,1] at high fractions or on noisy probes
  })
  colnames(purified_betas) <- sprintf("%s_pur%02d", tumor_sample_id, round(fractions * 100))
  rownames(purified_betas) <- common_probes
 
  # Project onto the SAME PCA space as before, restricted to common_probes -
  # this is a partial projection (subset of the original loadings), an
  # approximation rather than an exact reproduction of the full-probe
  # projection. Fine as long as common_probes covers most of probe_set
  # (see the warning above); degrades if coverage is low.
  rotation_sub <- pca_object$rotation[common_probes, , drop = FALSE]
  center_sub   <- pca_object$center[common_probes]
  scale_sub    <- if (isFALSE(pca_object$scale)) FALSE else pca_object$scale[common_probes]
  scores_new <- scale(t(purified_betas), center_sub, scale_sub) %*% rotation_sub
 
  data.frame(
    Sample_ID = rownames(scores_new),
    Source_Sample = tumor_sample_id,
    Cell_Type = cell_type,
    fraction_removed = fractions,
    n_probes_used = length(common_probes),
    as.data.frame(scores_new[, seq_len(min(5, ncol(scores_new))), drop = FALSE]),
    row.names = NULL
  )
}

ref_url <- "https://raw.githubusercontent.com/livjuergensen/InSilicoPurification/main/data/ref.csv"
ref_signatures <- read.csv(ref_url, row.names = 1, check.names = FALSE)
cat("Reference cell types available:", paste(colnames(ref_signatures), collapse = ", "), "\n")
cat(sprintf("Reference covers %d probes\n", nrow(ref_signatures)))


tumor_samples <- embeddings_PA %>% filter(family == 'Unclass.') %>% pull(Sample_ID) %>% unique()

cat(sprintf("Purifying %d unclassified samples\n", length(tumor_samples)))
ratios_to_use <- c(0,0.05, 0.1,0.15,0.2,0.25)  # matches your dilution series exactly;
                                            # 0 kept internally as the arrow start point
 
## ---- Build your own Sturm-derived reference signatures ---------------------
## Same three controls you already used for the additive dilution series -
## run through the identical purification pipeline as the five published
## signatures. Averaged across all samples matching each control type
## (matching the rowMeans() approach make_mixture_series() already uses for
## multiple normal_sample_ids), restricted to the PCA's own probe set so no
## extra probe-alignment logic is needed beyond what make_purification_series()
## already does internally.
sturm_reactive_samples    <- embeddings_all %>%
  filter(family == 'Control', Cohort == 'Sturm', PA == 'Non-neoplastic tissue',
         Methylation_Class == 'Control, reactive tissue') %>%
  pull(Sample_ID) %>% unique()
sturm_cerebellar_samples  <- embeddings_all %>%
  filter(family == 'Control', Cohort == 'Sturm', PA == 'Non-neoplastic tissue',
         Methylation_Class == 'Control, cerebellar tissue') %>%
  pull(Sample_ID) %>% unique()
sturm_hemispheric_samples <- embeddings_all %>%
  filter(family == 'Control', Cohort == 'Sturm', PA == 'Non-neoplastic tissue',
         Methylation_Class == 'Control, hemispheric tissue') %>%
  pull(Sample_ID) %>% unique()
 
cat(sprintf("Sturm controls found: reactive=%d, cerebellar=%d, hemispheric=%d\n",
            length(sturm_reactive_samples), length(sturm_cerebellar_samples),
            length(sturm_hemispheric_samples)))
stopifnot(length(sturm_reactive_samples) > 0, length(sturm_cerebellar_samples) > 0,
          length(sturm_hemispheric_samples) > 0)
 
pca_probe_set <- rownames(pca_full_PA$rotation)  # same probes the PCA was fit on
ref_signatures_sturm <- data.frame(
  Reactive    = rowMeans(beta_values_combat[pca_probe_set, sturm_reactive_samples, drop = FALSE]),
  Cerebellar  = rowMeans(beta_values_combat[pca_probe_set, sturm_cerebellar_samples, drop = FALSE]),
  Hemispheric = rowMeans(beta_values_combat[pca_probe_set, sturm_hemispheric_samples, drop = FALSE]),
  row.names = pca_probe_set,
  check.names = FALSE
)
 
purification_results_published <- do.call(rbind, lapply(colnames(ref_signatures), function(ct) {
  do.call(rbind, lapply(tumor_samples, function(s) {
    make_purification_series(
      tumor_sample_id = s, ref_signatures = ref_signatures, cell_type = ct,
      beta_matrix = beta_values_combat, pca_object = pca_full_PA,
      fractions = ratios_to_use
    )
  }))
})) %>% mutate(Reference_Source = "Published (Jurgensen et al.)")
 
purification_results_sturm <- do.call(rbind, lapply(colnames(ref_signatures_sturm), function(ct) {
  do.call(rbind, lapply(tumor_samples, function(s) {
    make_purification_series(
      tumor_sample_id = s, ref_signatures = ref_signatures_sturm, cell_type = ct,
      beta_matrix = beta_values_combat, pca_object = pca_full_PA,
      fractions = ratios_to_use
    )
  }))
})) %>% mutate(Reference_Source = "Sturm controls")
 
purification_results_matched <- rbind(purification_results_published, purification_results_sturm)
 
## ---- Facet labels, ordered, excluding the 0% start point -------------------
## Using "-X% normal" rather than "+X% normal" here deliberately, since this
## is the opposite direction from the additive dilution plot - removing
## reference signal, not adding it.
dilution_levels <- paste0('-', sprintf("%d%% normal", round(sort(unique(
  purification_results_matched$fraction_removed[purification_results_matched$fraction_removed > 0]
)) * 100)))
purification_results_matched <- purification_results_matched %>%
  mutate(dilution_label = factor(
    paste0('-', sprintf("%d%% normal", round(fraction_removed * 100))),
    levels = dilution_levels
  ))
 
## ---- Arrow starts: identical across cell types at fraction=0 --------------
arrow_starts_pur <- purification_results_matched %>%
  filter(fraction_removed == 0) %>%
  select(Source_Sample, PC1_start = PC1, PC2_start = PC2) %>%
  unique()   # defensive - same duplicate-row issue as before could recur
             # if embeddings_PA (or tumor_samples derived from it) still
             # has unresolved duplicates
 
arrow_df_pur <- purification_results_matched %>%
  filter(fraction_removed > 0) %>%
  left_join(arrow_starts_pur, by = "Source_Sample") %>%
  unique()
 
## ---- Combined plot, colored by reference cell type -------------------------
background_df <- embeddings_PA
all_cell_types <- c(colnames(ref_signatures), colnames(ref_signatures_sturm))
celltype_pal <-


    setNames(
  color_pal[c("Blue", "Orange", "Red", "Green", "Purple", "Teal", "Brown", "Pink")][seq_along(all_cell_types)],
  all_cell_types
)
print(celltype_pal)  # confirm the cell-type -> color mapping reads sensibly before trusting the plot
all_cell_types
celltype_pal <- color_pal[c('Red','Green','Orange','Teal')]
celltype_selection <-c('T cells','Microglia','Monocytes')
names(celltype_pal) <- celltype_selection

arrow_df_pur %>% head()

p_purification_matched <- ggplot() +
     geom_point(data = background_df %>% filter(family != 'PA'), aes(PC1, PC2),
             color = color_pal[['Dark Gray']], size = 2) +
    geom_point(data = background_df %>% filter(family == 'PA'), aes(PC1, PC2),
               color = color_pal[['Blue']], size = 2) +
  geom_segment(
    data = arrow_df_pur %>% filter(Cell_Type %in% celltype_selection,PC1_start<10 ),
    aes(x = PC1_start, y = PC2_start, xend = PC1, yend = PC2, color = Cell_Type,
        group = interaction(Source_Sample, Cell_Type)),
    arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
    alpha = 0.5, linewidth = 0.25
  ) +
  geom_point(data = arrow_df_pur %>% filter(Cell_Type %in% celltype_selection,PC1_start<10 ,fraction_removed != 0.05), aes(PC1, PC2, color = Cell_Type), size = 1, alpha = 0.5) +
  scale_color_manual(values = celltype_pal) +
  facet_wrap(~dilution_label, nrow=1) +
  theme_classic(base_size = 13) +
    labs(
         x = paste0('PC1 (',round(var_explained[1]),'%)' ),
        y = paste0('PC2 (',round(var_explained[2]),'%)' ), color = "Jürgensen et al reference") +
    theme(legend.position = 'bottom')
 
pdf('Purification_PA_published_and_sturm_references.pdf', height = 3, width = 8)
p_purification_matched
dev.off()





p_purification_M06 <-
    ggplot() +
     geom_point(data = background_df %>% filter(!grepl('M06', Sample_ID)), aes(PC1, PC2),
             color = color_pal[['Dark Gray']], size = 2) +
    geom_point(data = background_df %>% filter(grepl('M06', Sample_ID)), aes(PC1, PC2),
               color = color_pal[['Blue']], size = 2) +
  geom_segment(
    data = arrow_df_pur %>% filter(Cell_Type == 'Microglia',PC1_start<10 ,grepl('M06', Sample_ID)),
    aes(x = PC1_start, y = PC2_start, xend = PC1, yend = PC2, color = Cell_Type,
        group = interaction(Source_Sample, Cell_Type)),
    arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
    alpha = 0.5, linewidth = 0.25
  ) +
  geom_point(data = arrow_df_pur %>% filter(Cell_Type == 'Microglia',PC1_start<10,grepl('M06', Sample_ID)), aes(PC1, PC2, color = Cell_Type), size = 1, alpha = 0.5) +
  scale_color_manual(values = celltype_pal) +
  facet_wrap(~dilution_label, nrow=1) +
  theme_classic(base_size = 13) +
    labs(
         x = paste0('PC1 (',round(var_explained[1]),'%)' ),
        y = paste0('PC2 (',round(var_explained[2]),'%)' ), color = "Jürgensen et al reference") +
    theme(legend.position = 'bottom')



pdf('Purification_PA_published_and_Microglioa_MINT_M06only.pdf', height = 3, width = 8)
p_purification_M06
dev.off()





write.csv(purification_results_matched, "PA_purification_published_and_sturm.csv", row.names = FALSE)
 

arrow_df_pur$fraction_removed
















non_NF1 <- samplesheet_ref_pediatric %>% filter(!grepl('Neurofibromatosis', germline)) %>% pull(Sample_Name)
PA_samples_non_NF1_Sturm  <- embeddings_all %>% filter(Cohort == 'Sturm', Sample_ID %in% non_NF1) %>%
    filter(Subtype_Pathology == 'PA',family == 'PA') %>% pull(Sample_ID)

PA_samples_extended <- c(PA_samples,PA_samples_non_NF1_Sturm)

embeddings_PA_extended <- Create_embedding_df(mvalues_combat[,PA_samples_extended] , perplexity = 10) %>%
    mutate(Cohort = case_when(
               Cohort== 'Sturm' & Sample_ID %in% non_NF1 ~ 'Sturm NF1wt',
               Cohort == 'Sturm' & !Sample_ID %in% non_NF1 ~ 'Sturm NF1',
               TRUE ~ Cohort))
               


pca_full <- Run_PCA(mvalues_combat[,PA_samples_extended])
var_explained <- (pca_full$sdev^2) / sum(pca_full$sdev^2) * 100

pdf('PCA_PA_extended_gliomas_w_batch_correction_family.pdf', height = 3, width = 6)
PCA_PA_extended <- embeddings_PA_extended %>%
    arrange(desc(Cohort)) %>% 
    ggplot(aes(PC1, PC2, fill = family, shape = Cohort, color = Cohort)) +
    geom_point( size = 3,stroke = 1) +
    scale_fill_manual(values = subtype_pal) +
    scale_color_manual(values = cohort_pal) +
    scale_shape_manual(values = rev(c(22,24,21,23))) +
    theme_classic(base_size = 13) +
    labs(
        fill = 'DKFZ class',
        x = paste0('PC1 (',round(var_explained[1]),'%)' ),
        y = paste0('PC2 (',round(var_explained[2]),'%)' ))+
    guides(
        fill = guide_legend(override.aes = list(shape = 21))) +
    geom_point(data = embeddings_PA_extended %>% filter(Methylation_Class == 'PA, NF1-associated'),aes(PC1,PC2), shape = 8, color = 'black', inherit.aes = F, size = 1)
PCA_PA_extended
dev.off()
colnames(embeddings_PA_extended)


plot_data <- cbind(embeddings_PA_extended %>% select(Sample_ID,family,Cohort), pca_full$x[,2:20])

## ---- 3. Sample x sample Pearson correlation, based on PC scores ----------
pc_cols <- grep("^PC[0-9]+$", colnames(plot_data), value = TRUE) 
pc_matrix <- as.matrix(plot_data[, pc_cols])
rownames(pc_matrix) <- plot_data$Sample_ID
 
sample_cor <- cor(t(pc_matrix), method = "pearson")  # samples x samples
 
## ---- 4. Annotation data, in the SAME sample order as sample_cor ----------
annotation_df <- plot_data[match(rownames(sample_cor), plot_data$Sample_ID),
                            c("Cohort", "family")]


top_anno <- HeatmapAnnotation(
  Cohort = annotation_df$Cohort,
  `DKFZ class` = annotation_df$family,
  annotation_name_gp = gpar(fontsize = 8),
  col = list(Cohort = unlist(cohort_pal), `DKFZ class` = unlist(subtype_pal)),
  which = "column",
  show_legend = TRUE
)
left_anno <- HeatmapAnnotation(
  Cohort = annotation_df$Cohort,
  `DKFZ class` = annotation_df$family,
  annotation_name_gp = gpar(fontsize = 8),

  col = list(Cohort = unlist(cohort_pal), `DKFZ class` = unlist(subtype_pal)),
  which = "row",
  show_legend = FALSE  # avoid duplicating the same legend twice (top already shows it)
)
 
col_fun <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))
 
ht_sample_corr <- Heatmap(
  sample_cor,
  name = "Pearson r",
  col = col_fun,
  top_annotation = top_anno,
  left_annotation = left_anno,
  show_row_names = FALSE, show_column_names = FALSE,
  cluster_rows = TRUE, cluster_columns = TRUE,
  column_title = "Sample-sample correlation (PC2-PC50)"
)
 
## ---- 6. Draw standalone, to check it before combining ---------------------
pdf("Sample_correlation_heatmap.pdf", height = 8, width = 8)
draw(ht_sample_corr)
dev.off()
 

heatmap_grob <- grid.grabExpr(draw(ht_sample_corr, merge_legend = T))
combined_plot <- cowplot::plot_grid(PCA_PA_extended, heatmap_grob, labels = c('A','B'))
pdf("All_PA_PCA_and_sample_correlation_combined.pdf", height = 4.5, width = 10)
combined_plot
dev.off()



write.csv(embeddings_PA_extended %>%
            select(Sample_ID,Cohort,Platform,PC1,PC2,Loc,family,Subtype_Pathology),
            file = 'PA_extended_PCA_and_subtypes.csv',row.names = F, quote=F)


HGAP_samples <- embeddings_all_selected %>%
    filter(Cohort %in% c('MINT','Lucas','Sturm NF1'),
           family %in% c('HGAP','HGAP, NF1','Control') |
           (family == 'Unclass.' & Subtype_Pathology == 'HGAP') |
           (family == 'Control' & Subtype_Pathology == 'HGAP')) %>%
    pull(Sample_ID)

embeddings_HGAP <- Create_embedding_df(mvalues_combat[,HGAP_samples] , perplexity = 10) %>%
    mutate(Cohort = ifelse(Cohort == 'Sturm','Sturm NF1', Cohort)) %>%
    arrange(Cohort)
pca_full <- Run_PCA(mvalues_combat[,HGAP_samples])
var_explained <- (pca_full$sdev^2) / sum(pca_full$sdev^2) * 100

pdf('PCA_HGAP_gliomas_w_batch_correction_family.pdf', height = 3, width = 6)
PCA_HGAP <- embeddings_HGAP %>%
    ggplot(aes(PC1, PC2, fill = family, shape = Cohort, color = Cohort)) +
    geom_point( size = 4,stroke = 1) +
    scale_fill_manual(values = subtype_pal) +
    scale_color_manual(values = cohort_pal) +
    scale_shape_manual(values = c(22,24,21,23)) +
    theme_classic(base_size = 13) +
    labs(
        fill = 'DKFZ class',
        x = paste0('PC1 (',round(var_explained[1]),'%)' ),
        y = paste0('PC2 (',round(var_explained[2]),'%)' ))+
    guides(
        fill = guide_legend(override.aes = list(shape = 21))) +
    geom_point(data = embeddings_HGAP %>% filter(Methylation_Class == 'HGAP, NF1-associated'),aes(PC1,PC2), shape = 8, color = 'black', inherit.aes = F, size = 1)
PCA_HGAP
dev.off()

write.csv(embeddings_HGAP %>%
            select(Sample_ID,Cohort,Platform,PC1,PC2,Loc,family,Subtype_Pathology),
            file = 'HGAP_PCA_and_subtypes.csv',row.names = F, quote=F)


# Create
mixing_df <- data.frame(tumor = embeddings_HGAP %>% filter(family == 'HGAP') %>% pull(Sample_ID))
mixing_df$normal_reactive <- embeddings_all %>% filter(family == 'Control',Cohort == 'Sturm', PA == 'Non-neoplastic tissue', Methylation_Class == 'Control, reactive tissue') %>% pull(Sample_ID)
mixing_df$normal_cerebellar <- embeddings_all %>% filter(family == 'Control',Cohort == 'Sturm', PA == 'Non-neoplastic tissue', Methylation_Class == 'Control, cerebellar tissue') %>% pull(Sample_ID)
mixing_df$normal_hemispheric <- embeddings_all %>% filter(family == 'Control',Cohort == 'Sturm', PA == 'Non-neoplastic tissue', Methylation_Class == 'Control, hemispheric tissue') %>% pull(Sample_ID) %>% head(1)



mixing_results_all <- bind_rows(
  run_mixing_for_control("normal_reactive",    "Reactive"),
  run_mixing_for_control("normal_cerebellar",  "Cerebellar"),
  run_mixing_for_control("normal_hemispheric", "Hemispheric")
)
 
## ---- Facet labels, ordered, excluding the 0% start point -----------------
dilution_levels <- paste0('+', sprintf("%d%% normal", round(sort(unique(
  mixing_results_all$ratio_normal[mixing_results_all$ratio_normal > 0]
)) * 100)))
mixing_results_all <- mixing_results_all %>%
  mutate(dilution_label = factor(
    paste0('+', sprintf("%d%% normal", round(ratio_normal * 100))),
    levels = dilution_levels
  ))
 
## ---- Arrow starts: identical across control types at ratio=0 (the normal
## sample doesn't factor in yet), so one lookup per tumor sample suffices -
## matches what the original per-type scripts already did.
arrow_starts <- mixing_results_all %>%
  filter(ratio_normal == 0) %>%
  select(Source_Sample, PC1_start = PC1, PC2_start = PC2) %>%
  unique()
 
arrow_df_all <- mixing_results_all %>%
  filter(ratio_normal > 0) %>%
  left_join(arrow_starts, by = "Source_Sample") %>%
  unique()
 color_pal
## ---- Single combined plot, colored by control type ------------------------
background_df <- embeddings_HGAP
control_pal <- c(
  "Reactive"    = color_pal[['Red']],
  "Cerebellar"  = color_pal[['Green']],
  "Hemispheric" = color_pal[['Orange']]
)


p_dilution_combined <- ggplot() +
  geom_point(data = background_df %>% filter(family != 'HGAP'), aes(PC1, PC2),
             color = color_pal[['Dark Gray']], size = 2) +
    geom_point(data = background_df %>% filter(family == 'HGAP'), aes(PC1, PC2),
             color = color_pal[['Yellow']], size = 2) +

    
  geom_segment(
    data = arrow_df_all,
    aes(x = PC1_start, y = PC2_start, xend = PC1, yend = PC2, color = Control_Type,
        group = interaction(Source_Sample, Control_Type)),
    arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
    alpha = 0.8, linewidth = 0.6
  ) +
  geom_point(data = arrow_df_all, aes(PC1, PC2, color = Control_Type), size = 2, alpha = 0.6) +
  scale_color_manual(values = control_pal) +
  facet_wrap(~dilution_label) +
  theme_classic(base_size = 13) +
    labs(title = "In-silico normal-dilution by control tissue type",
         x = paste0('PC1 (',round(var_explained[1]),'%)' ),
        y = paste0('PC2 (',round(var_explained[2]),'%)' ), color = "Normal tissue")



pdf('Normal_dilutions_HGAP_combined.pdf', height = 6, width = 8)
p_dilution_combined
dev.off()








## =============================================================================
## In-silico PURIFICATION (subtraction), adapted from Jürgensen et al. 2025
## "In silico purification improves DNA methylation-based classification
## rates of pediatric low-grade gliomas" (Acta Neuropathol 150:34)
## github.com/livjuergensen/InSilicoPurification
##
## Mirror-image of make_mixture_series(): instead of ADDING a normal
## reference (diluting toward normal, which is what you did to explore the
## tumor->normal axis), this SUBTRACTS it, sweeping the assumed
## contamination fraction and re-projecting each purified profile onto your
## EXISTING PCA space - to see whether purification moves low-confidence
## (likely high-normal-contamination) samples back toward the high-
## confidence cluster, using the exact same linear demixing formula the
## paper's own method uses:
##
##   purified = (observed - f * reference) / (1 - f)
##
## swept across f (paper's default grid: 0 to 0.98, step 0.02).
## =============================================================================
 
## ---- Download the reference signatures ------------------------------------
## Median beta-values per non-malignant cell type (microglia, monocytes,
## neutrophils, T cells, neurons) from sorted/enriched cell populations -
## NOT verified against the live file structure here, so check the printed
## colnames() output below before trusting TARGET_CELL_TYPE further down.


## ---- Purification function (single cell type, single sample) -------------

embeddings_PA %>% filter(family == 'Unclass.' & PC1 < 10) %>% filter(Cohort == 'MINT')

samplesheet_MINT %>% filter(grepl('M06|M31|M32|M36',patient))
samplesheet %>% filter(grepl('M06|M31|M32|M36',patient))
tumor_samples <- embeddings_PA %>% filter(family == 'Unclass.' & PC1 < 10) %>% pull(Sample_ID) %>% unique()
ratios_to_use <- c(0,0.1,0.2,0.3,0.4 )  # matches your dilution series exactly;
                                            # 0 kept internally as the arrow start point
purification_results_matched <- do.call(rbind, lapply(colnames(ref_signatures), function(ct) {
  do.call(rbind, lapply(tumor_samples, function(s) {
    make_purification_series(
      tumor_sample_id = s, ref_signatures = ref_signatures, cell_type = ct,
      beta_matrix = beta_values_combat, pca_object = pca_full,
      fractions = ratios_to_use
    )
  }))
}))

## ---- Facet labels, ordered, excluding the 0% start point -------------------
## Using "-X% normal" rather than "+X% normal" here deliberately, since this
## is the opposite direction from the additive dilution plot - removing
## reference signal, not adding it.
dilution_levels <- paste0('-', sprintf("%d%% normal", round(sort(unique(
  purification_results_matched$fraction_removed[purification_results_matched$fraction_removed > 0]
)) * 100)))
purification_results_matched <- purification_results_matched %>%
  mutate(dilution_label = factor(
    paste0('-', sprintf("%d%% normal", round(fraction_removed * 100))),
    levels = dilution_levels
  ))
 
## ---- Arrow starts: identical across cell types at fraction=0 --------------
arrow_starts_pur <- purification_results_matched %>%
  filter(fraction_removed == 0) %>%
  select(Source_Sample, PC1_start = PC1, PC2_start = PC2) %>%
  unique()   # defensive - same duplicate-row issue as before could recur
             # if embeddings_PA (or tumor_samples derived from it) still
             # has unresolved duplicates
 
arrow_df_pur <- purification_results_matched %>%
  filter(fraction_removed > 0) %>%
  left_join(arrow_starts_pur, by = "Source_Sample") %>%
  unique()
 
## ---- Combined plot, colored by reference cell type -------------------------
background_df <- embeddings_PA
celltype_pal <- setNames(
  color_pal[c( "Orange", "Red", "Green", "Purple",'Teal')][seq_along(colnames(ref_signatures))],
  colnames(ref_signatures)
)



p_purification_matched <- ggplot() +
     geom_point(data = background_df %>% filter(family != 'PA'), aes(PC1, PC2),
             color = color_pal[['Dark Gray']], size = 2) +
    geom_point(data = background_df %>% filter(family == 'PA'), aes(PC1, PC2),
               color = color_pal[['Blue']], size = 2) +
  geom_segment(
    data = arrow_df_pur,
    aes(x = PC1_start, y = PC2_start, xend = PC1, yend = PC2, color = Cell_Type,
        group = interaction(Source_Sample, Cell_Type)),
    arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
    alpha = 0.75, linewidth = 0.5
  ) +
  geom_point(data = arrow_df_pur, aes(PC1, PC2, color = Cell_Type), size = 1.5, alpha = 0.5) +
  scale_color_manual(values = celltype_pal) +
  facet_wrap(~dilution_label) +
  theme_classic(base_size = 13) +
    labs(
         x = paste0('PC1 (',round(var_explained[1]),'%)' ),
        y = paste0('PC2 (',round(var_explained[2]),'%)' ), color = "Reference cell type")
 
pdf('Purification_PA_reference_celltypes_matched_ratios.pdf', height = 6, width = 8)
p_purification_matched
dev.off()
 
write.csv(purification_results_matched, "PA_purification_matched_ratios.csv", row.names = FALSE)


tumor_samples <- embeddings_PA_extended %>% filter(family == 'Unclass.') %>% pull(Sample_ID) %>% unique()
ratios_to_use <- c(0,0.1,0.2,0.3,0.4 )  # matches your dilution series exactly;
                                            # 0 kept internally as the arrow start point
purification_results_matched <- do.call(rbind, lapply(colnames(ref_signatures), function(ct) {
  do.call(rbind, lapply(tumor_samples, function(s) {
    make_purification_series(
      tumor_sample_id = s, ref_signatures = ref_signatures, cell_type = ct,
      beta_matrix = beta_values_combat, pca_object = pca_full,
      fractions = ratios_to_use
    )
  }))
}))

## ---- Facet labels, ordered, excluding the 0% start point -------------------
## Using "-X% normal" rather than "+X% normal" here deliberately, since this
## is the opposite direction from the additive dilution plot - removing
## reference signal, not adding it.
dilution_levels <- paste0('-', sprintf("%d%% normal", round(sort(unique(
  purification_results_matched$fraction_removed[purification_results_matched$fraction_removed > 0]
)) * 100)))
purification_results_matched <- purification_results_matched %>%
  mutate(dilution_label = factor(
    paste0('-', sprintf("%d%% normal", round(fraction_removed * 100))),
    levels = dilution_levels
  ))
 
## ---- Arrow starts: identical across cell types at fraction=0 --------------
arrow_starts_pur <- purification_results_matched %>%
  filter(fraction_removed == 0) %>%
  select(Source_Sample, PC1_start = PC1, PC2_start = PC2) %>%
  unique()   # defensive - same duplicate-row issue as before could recur
             # if embeddings_PA_extended (or tumor_samples derived from it) still
             # has unresolved duplicates

arrow_df_pur <- purification_results_matched %>%
  filter(fraction_removed > 0) %>%
  left_join(arrow_starts_pur, by = "Source_Sample") %>%
  unique()
 
## ---- Combined plot, colored by reference cell type -------------------------
background_df <- embeddings_PA_extended
celltype_pal <- setNames(
  color_pal[c( "Orange", "Red", "Green", "Purple",'Teal')][seq_along(colnames(ref_signatures))],
  colnames(ref_signatures)
)

background_df %>% filter(family == 'Unclass.')


p_purification_matched <- ggplot() +
     geom_point(data = background_df %>% filter(family != 'PA'), aes(PC1, PC2),
             color = color_pal[['Dark Gray']], size = 2) +
    geom_point(data = background_df %>% filter(family == 'PA'), aes(PC1, PC2),
               color = color_pal[['Blue']], size = 2) +
  geom_segment(
    data = arrow_df_pur,
    aes(x = PC1_start, y = PC2_start, xend = PC1, yend = PC2, color = Cell_Type,
        group = interaction(Source_Sample, Cell_Type)),
    arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
    alpha = 0.75, linewidth = 0.5
  ) +
  geom_point(data = arrow_df_pur, aes(PC1, PC2, color = Cell_Type), size = 1.5, alpha = 0.5) +
  scale_color_manual(values = celltype_pal) +
  facet_wrap(~dilution_label) +
  theme_classic(base_size = 13) +
  labs(title = "In-silico purification using Jürgensen et al. reference signatures",
       x = "PC1", y = "PC2", color = "Reference cell type")
 
pdf('Purification_PA_extended_reference_celltypes_matched_ratios.pdf', height = 6, width = 8)
p_purification_matched
 dev.off()




tumor_samples <- embeddings_HGAP %>% filter(family == 'Unclass.') %>% pull(Sample_ID) %>% unique()
ratios_to_use <- c(0,0.1,0.2,0.3,0.4 )  # matches your dilution series exactly;
                                            # 0 kept internally as the arrow start point
purification_results_matched <- do.call(rbind, lapply(colnames(ref_signatures), function(ct) {
  do.call(rbind, lapply(tumor_samples, function(s) {
    make_purification_series(
      tumor_sample_id = s, ref_signatures = ref_signatures, cell_type = ct,
      beta_matrix = beta_values_combat, pca_object = pca_full,
      fractions = ratios_to_use
    )
  }))
}))

## ---- Facet labels, ordered, excluding the 0% start point -------------------
## Using "-X% normal" rather than "+X% normal" here deliberately, since this
## is the opposite direction from the additive dilution plot - removing
## reference signal, not adding it.
dilution_levels <- paste0('-', sprintf("%d%% normal", round(sort(unique(
  purification_results_matched$fraction_removed[purification_results_matched$fraction_removed > 0]
)) * 100)))
purification_results_matched <- purification_results_matched %>%
  mutate(dilution_label = factor(
    paste0('-', sprintf("%d%% normal", round(fraction_removed * 100))),
    levels = dilution_levels
  ))
 
## ---- Arrow starts: identical across cell types at fraction=0 --------------
arrow_starts_pur <- purification_results_matched %>%
  filter(fraction_removed == 0) %>%
  select(Source_Sample, PC1_start = PC1, PC2_start = PC2) %>%
  unique()   # defensive - same duplicate-row issue as before could recur
             # if embeddings_HGAP (or tumor_samples derived from it) still
             # has unresolved duplicates

arrow_df_pur <- purification_results_matched %>%
  filter(fraction_removed > 0) %>%
  left_join(arrow_starts_pur, by = "Source_Sample") %>%
  unique()
 
## ---- Combined plot, colored by reference cell type -------------------------
background_df <- embeddings_HGAP
celltype_pal <- setNames(
  color_pal[c( "Orange", "Red", "Green", "Purple",'Teal')][seq_along(colnames(ref_signatures))],
  colnames(ref_signatures)
)

background_df %>% filter(family == 'Unclass.')


p_purification_matched <- ggplot() +
     geom_point(data = background_df %>% filter(family != 'HGAP'), aes(PC1, PC2),
             color = color_pal[['Dark Gray']], size = 2) +
    geom_point(data = background_df %>% filter(family == 'HGAP'), aes(PC1, PC2),
               color = color_pal[['Yellow']], size = 2) +
  geom_segment(
    data = arrow_df_pur,
    aes(x = PC1_start, y = PC2_start, xend = PC1, yend = PC2, color = Cell_Type,
        group = interaction(Source_Sample, Cell_Type)),
    arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
    alpha = 0.75, linewidth = 0.5
  ) +
  geom_point(data = arrow_df_pur, aes(PC1, PC2, color = Cell_Type), size = 1.5, alpha = 0.5) +
  scale_color_manual(values = celltype_pal) +
  facet_wrap(~dilution_label) +
  theme_classic(base_size = 13) +
  labs(
       x = "PC1", y = "PC2", color = "Jürgensen et al reference")
 
pdf('Purification_HGAP_reference_celltypes_matched_ratios.pdf', height = 6, width = 8)
p_purification_matched
 dev.off()
 
write.csv(purification_results_matched, "HGAP_purification_matched_ratios.csv", row.names = FALSE)

















# =============================================================================
## LUMP (Leukocytes UnMethylation for Purity) - Aran, Sirota & Butte 2015
## Nat Commun 6:8971. Formula: purity = mean(beta at 44 CpGs) / 0.85
##
## These 44 CpGs are consistently UNmethylated (<5%) in leukocytes and
## consistently methylated (>30%) across 21 TCGA cancer types - so a HIGH
## average beta at these sites means mostly non-immune (tumor) cells; a LOW
## average means mostly immune infiltration. The formula gives PURITY
## directly, not immune fraction - immune % is the complement (1 - purity).
##
## NOTE ON THE PROBE LIST BELOW: your pasted IDs had inconsistent digit
## counts (cg240653 = 6 digits, cg1138020 = 7 digits, cg10511890 = 8 digits) -
## real Illumina probe IDs are always exactly 8 digits after "cg". This is
## the classic leading-zeros-stripped-by-Excel/PDF-extraction corruption.
## Fixed by zero-padding every numeric suffix to 8 digits below. Two
## non-probe fragments from figure-caption text ("Consistently unmethylated
## sites...", "Consistently methylated sites...") were also excluded - the
## remaining count is exactly 44, matching the paper.
## =============================================================================
 
lump_cpgs_raw <- c(
  "cg240653","cg450164","cg880290","cg933696","cg1138020","cg2026204",
  "cg2053964","cg2167021","cg2997560","cg3431741","cg3436397","cg3841065",
  "cg4915566","cg5199874","cg5305434","cg5769344","cg5798125","cg7002058",
  "cg7598052","cg7641284","cg8854008","cg9302355","cg9606470","cg10511890",
  "cg10559416","cg13030790","cg13912307","cg14076977","cg14913777","cg17518965",
  "cg19466818","cg20170223","cg20695297","cg21164509","cg21376733","cg22331159",
  "cg23114964","cg23553480","cg24796554","cg25384897","cg25574765","cg26427109",
  "cg26842802","cg27215100"
)
stopifnot(length(lump_cpgs_raw) == 44)
 
# Zero-pad every numeric suffix to 8 digits (standard Illumina cg-ID format)
lump_cpgs <- sprintf("cg%08d", as.integer(sub("^cg", "", lump_cpgs_raw)))
cat("Corrected LUMP probe IDs (first 5 shown):\n")
print(head(lump_cpgs, 5))
 
## ---- LUMP scoring function -------------------------------------------------
lump_score <- function(beta_matrix, cpg_ids = lump_cpgs) {
  present <- intersect(cpg_ids, rownames(beta_matrix))
  cat(sprintf("[LUMP] Using %d/%d published CpGs present in your data\n",
              length(present), length(cpg_ids)))
  if (length(present) == 0) {
    stop("None of the 44 LUMP CpGs were found in beta_matrix - check probe ",
         "ID formatting or whether these CpGs survived your QC/blacklist filtering.")
  }
  if (length(present) < length(cpg_ids) * 0.5) {
    warning("Fewer than half of the LUMP CpGs are present - estimates may be ",
            "less reliable than the original 44-site calculation.")
  }
 
  mean_meth <- colMeans(beta_matrix[present, , drop = FALSE], na.rm = TRUE)
  purity <- pmin(mean_meth / 0.85, 1)  # clip at 1 - the /0.85 calibration
                                        # constant can occasionally push a
                                        # very heavily-methylated sample
                                        # slightly over 1 otherwise
 
  data.frame(
    Sample_ID = names(mean_meth),
    LUMP_Tumor_pct = purity * 100,
    LUMP_Immune_pct = (1 - purity) * 100,
    n_cpgs_used = length(present),
    row.names = NULL
  )
}


## ---- Apply to your PA samples ---------------------------------------------
pa_betas <- beta_values_combat[, PA_samples, drop = FALSE]  # swap in whichever final beta matrix you want

pa_betas %>% head()
write.csv(lump_results, "PA_LUMP_purity_immune_pct.csv", row.names = FALSE)




## ---- Sanity-check plot: distribution by cohort -----------------------------
pdf('PCA1vs_LUMP_PA_gliomas_w_batch_correction_family.pdf', height = 3, width = 6)
PCA_Lump <- embeddings_PA %>%
    left_join(lump_results) %>%
    ggplot(aes(PC1, LUMP_Immune_pct, fill = family, shape = Cohort, color = Cohort)) +
    geom_point( size = 4,stroke = 1) +
    scale_fill_manual(values = subtype_pal) +
    scale_color_manual(values = cohort_pal) +
    scale_shape_manual(values = c(22,24,21,23)) +
    theme_classic(base_size = 13) +
    labs(
        fill = 'DKFZ class',
        x = paste0('PC1 (',round(var_explained[1]),'%)' ),
        y = 'LUMP leukocyte %-age')+
    guides(
        fill = guide_legend(override.aes = list(shape = 21))) 
PCA_Lump
dev.off()


pdf('LUMP_plot.pdf', height = 6, width = 6)
cowplot::plot_grid(PCA_PA,PCA_Lump, nrow = 2)
dev.off()


embeddings_PA <- embeddings_PA %>% left_join(lump_results)


## =============================================================================
## PC1/PC2 correlation dot-heatmaps (ComplexHeatmap)
## Color = Pearson correlation coefficient, Size = significance (-log10 p
## from cor.test)
##
## Two separate analyses, same visual encoding:
##   1. SAMPLE level:  PC1/PC2 scores  vs  sample variables (Age, Horvath
##      clock, etc.) - correlated ACROSS SAMPLES
##   2. LOADINGS level: PC1/PC2 loadings (per-probe weights)  vs  reference
##      cell-type methylation profiles (published + Sturm) - correlated
##      ACROSS PROBES
## =============================================================================

library(ComplexHeatmap)
library(circlize)
library(grid)

## =============================================================================
## Generic dot-heatmap helper - reused for both analyses below
## =============================================================================
plot_corr_dotplot <- function(cor_mat, pval_mat, title = "",
                               max_neglogp_cap = 10,
                               p_adjust_method = "BH",
                               use_adjusted_p = FALSE,
                               padj_border_threshold = 0.05,
                               border_color_sig = "red",
                               border_color_nonsig = "black",
                               dot_mm_scale = 6,
                               draw_now = TRUE) {
  stopifnot(identical(dim(cor_mat), dim(pval_mat)))
 
  col_fun <- colorRamp2(c(-1, 0, 1), c("#2166AC", "white", "#B2182B"))
 
  # Using RAW p-values by default now (use_adjusted_p = FALSE) - set to TRUE
  # to switch back to BH/Bonferroni-corrected values. Whichever is chosen
  # becomes the SOLE significance measure driving both dot size and border.
  padj_mat <- if (use_adjusted_p) {
    matrix(
      p.adjust(as.vector(pval_mat), method = p_adjust_method),
      nrow = nrow(pval_mat), ncol = ncol(pval_mat), dimnames = dimnames(pval_mat)
    )
  } else {
    pval_mat
  }
 
  # Significance -> dot size, via -log10(p), capped
  neglogp <- -log10(padj_mat)
  neglogp[!is.finite(neglogp)] <- max_neglogp_cap
  neglogp <- pmin(neglogp, max_neglogp_cap)
 
  min_r <- 0.1; max_r <- 0.48
  # FIXED absolute domain (0 to max_neglogp_cap), not data-relative
  # (range(neglogp) within this specific dataset). This matters once the
  # legend uses fixed breakpoints (p = 1, 0.01, 1e-5, 1e-10 below): with a
  # data-relative domain, a legend entry at a p-value that happens to fall
  # outside THIS dataset's observed range would extrapolate to a radius
  # that doesn't match any real dot in the plot - same class of mismatch
  # as the mm-scaling issue fixed earlier, just in the size mapping instead
  # of the unit system. A fixed domain guarantees any given p-value always
  # maps to the same dot size, in the plot or the legend, regardless of
  # what else is in the dataset.
  radius_mat <- if (max_neglogp_cap == 0) {
    matrix(mean(c(min_r, max_r)), nrow(neglogp), ncol(neglogp), dimnames = dimnames(neglogp))
  } else {
    min_r + (neglogp / max_neglogp_cap) * (max_r - min_r)
  }
 
  # Border color: black by default, red where p is below the threshold
  border_mat <- matrix(
    ifelse(!is.na(padj_mat) & padj_mat < padj_border_threshold,
           border_color_sig, border_color_nonsig),
    nrow = nrow(padj_mat), ncol = ncol(padj_mat), dimnames = dimnames(padj_mat)
  )
 
  # IMPORTANT: dot_mm_scale is used for BOTH the heatmap cells below AND the
  # legend further down - this is what guarantees a dot in the plot and a
  # dot in the legend at the same radius_mat value render at the SAME
  # physical size. Previously the heatmap used a cell-relative size
  # (radius * min(cell width, cell height)) while the legend used an
  # independent absolute-mm size - two different coordinate systems that
  # were never actually calibrated to agree, just coincidentally looked
  # similar under the old adjusted-p distribution. Switching BOTH to the
  # same absolute-mm formula removes that whole class of mismatch by
  # construction, rather than trying to reverse-engineer what the
  # cell-relative size would be in mm (which depends on final rendered
  # figure dimensions, not knowable cleanly at legend-construction time).
  ht <- Heatmap(
    cor_mat, name = "Pearson r",
    col = col_fun,
    rect_gp = gpar(type = "none"),
    cell_fun = function(j, i, x, y, w, h, fill) {
      grid.rect(x = x, y = y, width = w, height = h,
                gp = gpar(col = "grey90", fill = NA))
      if (!is.na(cor_mat[i, j])) {
        grid.circle(
          x = x, y = y, r = unit(radius_mat[i, j] * dot_mm_scale, "mm"),
          gp = gpar(fill = col_fun(cor_mat[i, j]), col = border_mat[i, j], lwd = 1.2)
        )
      }
    },
    cluster_rows = TRUE, cluster_columns = FALSE,
    show_row_names = TRUE, show_column_names = TRUE,
    column_title = title,
    row_names_gp = gpar(fontsize = 8),
    row_names_side = "left",
    heatmap_legend_param = list(title = "Pearson r", at = c(-1, -0.5, 0, 0.5, 1))
  )
 
  # Legend dots now use the IDENTICAL formula (radius_mat * dot_mm_scale)
  # as the heatmap cells above - guaranteed match, not a separate estimate.
  # Fixed breakpoints at p = 1, 0.01, 1e-5, 1e-10 (i.e. -log10(p) = 0, 2, 5, 10)
  # rather than automatically-generated pretty() breaks - fewer, cleaner
  # legend entries at round p-value landmarks.
  legend_breaks <- c(0, 2, 5, 10)
  legend_breaks <- legend_breaks[legend_breaks <= max_neglogp_cap]
  legend_radii <- if (max_neglogp_cap == 0) {
    rep(mean(c(min_r, max_r)), length(legend_breaks))
  } else {
    min_r + (legend_breaks / max_neglogp_cap) * (max_r - min_r)
  }
  p_label <- if (use_adjusted_p) sprintf("padj (%s)", p_adjust_method) else "p"
  size_legend <- Legend(
    labels = sprintf("%s = %s", p_label,
                      ifelse(legend_breaks == 0, "1", format(10^-legend_breaks, digits = 2, scientific = TRUE))),
    title = "Significance",
    type = "points", pch = 16,
    legend_gp = gpar(col = "black"),
    size = unit(legend_radii * dot_mm_scale, "mm"),
    grid_height = unit(max_r * dot_mm_scale + 2, "mm"),
    grid_width = unit(max_r * dot_mm_scale + 2, "mm")
  )
 
  border_legend <- Legend(
    labels = c(sprintf("%s < %s", p_label, format(padj_border_threshold, scientific = TRUE)),
               sprintf("%s >= %s", p_label, format(padj_border_threshold, scientific = TRUE))),
    title = "Dot border",
    type = "points", pch = 21,
    legend_gp = gpar(col = c(border_color_sig, border_color_nonsig), fill = "white"),
    size = unit(3, "mm")
  )
 
  legend_list <- list(size_legend, border_legend)
 
  # draw_now=TRUE (default) preserves the original behavior: draws
  # immediately, useful for simple standalone pdf(); plot_corr_dotplot();
  # dev.off() usage. Set draw_now=FALSE to get the constructed objects back
  # WITHOUT drawing - required for grid.grabExpr()/cowplot composition,
  # since grid.grabExpr() only captures drawing that happens INSIDE the
  # expression you pass it, not a value that was already drawn earlier as
  # a side effect of calling this function.
  if (draw_now) {
    draw(ht, annotation_legend_list = legend_list, merge_legend = TRUE)
  } else {
    list(ht = ht, legend_list = legend_list)
  }
}
 
## =============================================================================
## Unified Pearson-correlation function: continuous variables directly,
## categorical variables one-hot encoded
## =============================================================================
## Returns one row (or several, for multi-level categoricals) per variable,
## each with a Pearson r and a p-value from cor.test() - the SAME statistic
## throughout, never a different method depending on variable type:
##   - numeric target -> Pearson correlation directly
##   - factor/character, 2 levels -> Pearson correlation on the 0/1-coded
##     indicator (point-biserial correlation - mathematically identical to
##     Pearson r, just a name for the special case of one binary variable)
##     A binary variable produces ONE row (two levels would just be mirror
##     images of each other, r and -r with the same p-value - redundant).
##     A multi-level variable produces ONE ROW PER LEVEL (one-vs-rest: that
##     level coded 1, everything else coded 0), since a single omnibus
##     number can't say which specific level drives an association.
compute_pc_pearson_rows <- function(pc_vector, var_vector, var_name) {
  if (is.numeric(var_vector)) {
    complete_idx <- !is.na(pc_vector) & !is.na(var_vector)
    if (sum(complete_idx) < 3) {
      return(data.frame(row_label = var_name, estimate = NA_real_, p.value = NA_real_))
    }
    ct <- cor.test(pc_vector[complete_idx], var_vector[complete_idx], method = "pearson")
    return(data.frame(row_label = var_name, estimate = unname(ct$estimate), p.value = ct$p.value))
  }
 
  var_factor <- droplevels(as.factor(var_vector))
  n_levels <- nlevels(var_factor)
 
  if (n_levels < 2) {
    return(data.frame(row_label = var_name, estimate = NA_real_, p.value = NA_real_))
  }
 
  if (n_levels == 2) {
    indicator <- as.numeric(var_factor) - 1
    complete_idx <- !is.na(pc_vector) & !is.na(indicator)
    if (sum(complete_idx) < 3) {
      return(data.frame(row_label = var_name, estimate = NA_real_, p.value = NA_real_))
    }
    ct <- cor.test(pc_vector[complete_idx], indicator[complete_idx], method = "pearson")
    return(data.frame(
      row_label = sprintf("%s (%s vs %s)", var_name, levels(var_factor)[2], levels(var_factor)[1]),
      estimate = unname(ct$estimate), p.value = ct$p.value
    ))
  }
 
  do.call(rbind, lapply(levels(var_factor), function(lvl) {
    indicator <- as.numeric(var_factor == lvl)
    complete_idx <- !is.na(pc_vector) & !is.na(indicator)
    if (sum(complete_idx) < 3 || length(unique(indicator[complete_idx])) < 2) {
      return(data.frame(row_label = sprintf("%s = %s", var_name, lvl),
                         estimate = NA_real_, p.value = NA_real_))
    }
    ct <- cor.test(pc_vector[complete_idx], indicator[complete_idx], method = "pearson")
    data.frame(row_label = sprintf("%s = %s", var_name, lvl),
               estimate = unname(ct$estimate), p.value = ct$p.value)
  }))
}
 
## =============================================================================
## Unified Pearson-correlation function: continuous variables directly,
## categorical variables one-hot encoded
## =============================================================================
## Returns one row (or several, for multi-level categoricals) per variable,
## each with a Pearson r and a p-value from cor.test() - the SAME statistic
## throughout, never a different method depending on variable type:
##   - numeric target -> Pearson correlation directly
##   - factor/character, 2 levels -> Pearson correlation on the 0/1-coded
##     indicator (point-biserial correlation - mathematically identical to
##     Pearson r, just a name for the special case of one binary variable)
##     Pearson correlation, applied uniformly - continuous variables
##     correlated directly; categorical variables one-hot encoded first.
##     A binary variable produces ONE row (two levels would just be mirror
##     images of each other, r and -r with the same p-value - redundant).
##     A multi-level variable produces ONE ROW PER LEVEL (one-vs-rest: that
##     level coded 1, everything else coded 0), since a single omnibus
##     number can't say which specific level drives an association -
##     point-biserial correlation IS just Pearson correlation on the 0/1
##     indicator, so this is genuinely the same statistic throughout, not
##     three different methods dressed up to look similar. Whole approach
##     is one sentence to describe: "Pearson correlation; categorical
##     variables were one-hot encoded."
compute_pc_pearson_rows <- function(pc_vector, var_vector, var_name) {
  if (is.numeric(var_vector)) {
    complete_idx <- !is.na(pc_vector) & !is.na(var_vector)
    if (sum(complete_idx) < 3) {
      return(data.frame(row_label = var_name, estimate = NA_real_, p.value = NA_real_))
    }
    ct <- cor.test(pc_vector[complete_idx], var_vector[complete_idx], method = "pearson")
    return(data.frame(row_label = var_name, estimate = unname(ct$estimate), p.value = ct$p.value))
  }
 
  var_factor <- droplevels(as.factor(var_vector))
  n_levels <- nlevels(var_factor)
 
  if (n_levels < 2) {
    return(data.frame(row_label = var_name, estimate = NA_real_, p.value = NA_real_))
  }
 
  if (n_levels == 2) {
    # One row only - the second level vs the first (alphabetical by
    # default R factor ordering). Two mirrored rows would add nothing.
    indicator <- as.numeric(var_factor) - 1
    complete_idx <- !is.na(pc_vector) & !is.na(indicator)
    if (sum(complete_idx) < 3) {
      return(data.frame(row_label = var_name, estimate = NA_real_, p.value = NA_real_))
    }
    ct <- cor.test(pc_vector[complete_idx], indicator[complete_idx], method = "pearson")
    return(data.frame(
      row_label = sprintf("%s (%s vs %s)", var_name, levels(var_factor)[2], levels(var_factor)[1]),
      estimate = unname(ct$estimate), p.value = ct$p.value
    ))
  }
 
  # >2 levels: one-hot encode, one row per level (one-vs-rest)
  do.call(rbind, lapply(levels(var_factor), function(lvl) {
    indicator <- as.numeric(var_factor == lvl)
    complete_idx <- !is.na(pc_vector) & !is.na(indicator)
    if (sum(complete_idx) < 3 || length(unique(indicator[complete_idx])) < 2) {
      return(data.frame(row_label = sprintf("%s = %s", var_name, lvl),
                         estimate = NA_real_, p.value = NA_real_))
    }
    ct <- cor.test(pc_vector[complete_idx], indicator[complete_idx], method =cell "pearson")
    data.frame(row_label = sprintf("%s = %s", var_name, lvl),
               estimate = unname(ct$estimate), p.value = ct$p.value)
  }))
}
 

## =============================================================================
## 1. SAMPLE-LEVEL: PC1/PC2 vs Age, Horvath clock, Cohort, family, etc.
## =============================================================================
## Assumes Age and Horvath_Age columns already exist on embeddings_PA (or
## whichever data.frame you're using) - see Horvath_clock_calculation.R if
## not. Every variable in sample_vars is handled automatically by
## compute_pc_pearson_rows() above, regardless of type - no need to
## pre-process, dummy-code, or separate continuous from categorical yourself.
sample_vars <- c( "Horvath", "Cohort", "Platform", "LUMP_Immune_pct", 'Loc','Mean_methylation')
                 # extend freely - a multi-level variable expands into
                 # multiple rows automatically (one per level); add
                 # NF1_status, Subtype_Pathology, LUMP_Tumor_pct, etc. as
                 # relevant
pc_vars <- c("PC1", "PC2")

missing_cols <- setdiff(c(sample_vars, pc_vars), colnames(embeddings_PA))
if (length(missing_cols) > 0) {
  stop("Missing column(s) in embeddings_PA: ", paste(missing_cols, collapse = ", "),
       " - add these before running the sample-level heatmap.")
}

sample_df <- embeddings_PA %>% select(all_of(c(pc_vars, sample_vars)))

sample_results <- do.call(rbind, lapply(pc_vars, function(pv) {
  do.call(rbind, lapply(sample_vars, function(sv) {
    res <- compute_pc_pearson_rows(sample_df[[pv]], sample_df[[sv]], sv)
    res$pc <- pv
    res
  }))
}))

## Reshape into the cor_mat/pval_mat shape plot_corr_dotplot() expects.
## Row order follows sample_vars/level order as first encountered, so
## multi-level variables' rows stay grouped together rather than scattering.
row_labels <- unique(sample_results$row_label)

cor_mat_sample  <- matrix(NA_real_, nrow = length(row_labels), ncol = length(pc_vars),
                           dimnames = list(row_labels, pc_vars))
pval_mat_sample <- cor_mat_sample
for (k in seq_len(nrow(sample_results))) {
  r <- sample_results[k, ]
  cor_mat_sample[r$row_label, r$pc]  <- r$estimate
  pval_mat_sample[r$row_label, r$pc] <- r$p.value
}

rownames(cor_mat_sample)[rownames(cor_mat_sample) == 'LUMP_Immune_pct'] <- 'LUMP leukocyte %-age'
rownames(cor_mat_sample)[rownames(cor_mat_sample) == 'Mean_methylation'] <- 'Average m-value'
    

pdf("PC_correlation_dotheatmap_sample_level.pdf", height = 3.5, width = 5)
plot_corr_dotplot(cor_mat_sample, pval_mat_sample)
dev.off()

## =============================================================================
## 2. LOADINGS-LEVEL: PC1/PC2 loadings vs reference cell-type signatures
## =============================================================================
## Correlates each PC's per-probe LOADING (its weight/direction in pca_full$rotation)
## against each reference's per-probe methylation level, across probes -
## i.e. "does the PC1 axis look like this cell type's methylation pattern?"
## Uses both the published (Jurgensen et al.) and Sturm-derived references
## you already built.
loading_probes <- rownames(pca_full_PA$rotation)

ref_combined <- cbind(
  ref_signatures[loading_probes, , drop = FALSE],
  ref_signatures_sturm[loading_probes, , drop = FALSE]
)


rownames(ref_combined) <- loading_probes
ref_vars <- colnames(ref_combined)

cor_mat_loadings  <- matrix(NA_real_, nrow = length(ref_vars), ncol = length(pc_vars),
                             dimnames = list(ref_vars, pc_vars))
pval_mat_loadings <- cor_mat_loadings

for (pv in pc_vars) {
  loading_vec <- pca_full_PA$rotation[loading_probes, pv]
  for (rv in ref_vars) {
    ref_vec <- ref_combined[, rv]
    valid <- !is.na(ref_vec) & !is.na(loading_vec)
    if (sum(valid) < 3) {
      warning(sprintf("Fewer than 3 valid probes for %s vs %s - leaving as NA", rv, pv))
      next
    }
    ct <- cor.test(loading_vec[valid], ref_vec[valid], method = "pearson")
    cor_mat_loadings[rv, pv]  <- unname(ct$estimate)
    pval_mat_loadings[rv, pv] <- ct$p.value
  }
}
p_corr_loading <-   draw(plot_corr_dotplot(cor_mat_loadings, pval_mat_loadings ))
dev.off()
pdf("PC_correlation_dotheatmap_loadings_level.pdf", height = 3.5, width = 5)
p_corr_loading
dev.off()



## ---- Save the underlying numbers too, not just the figures -----------------
write.csv(as.data.frame(cor_mat_sample) %>% tibble::rownames_to_column("Variable"),
          "PC_correlation_sample_level_r.csv", row.names = FALSE)
write.csv(as.data.frame(pval_mat_sample) %>% tibble::rownames_to_column("Variable"),
          "PC_correlation_sample_level_pval.csv", row.names = FALSE)
write.csv(as.data.frame(cor_mat_loadings) %>% tibble::rownames_to_column("Reference"),
          "PC_correlation_loadings_level_r.csv", row.names = FALSE)
write.csv(as.data.frame(pval_mat_loadings) %>% tibble::rownames_to_column("Reference"),
          "PC_correlation_loadings_level_pval.csv", row.names = FALSE)






cowplot::plot_grid(PCA_PA,plot_corr_dotplot(cor_mat_sample, pval_mat_sample,draw_now = F))
class(heatmap_grob)

built <- plot_corr_dotplot(cor_mat_sample, pval_mat_sample, draw_now = FALSE)
heatmap_grob <- grid.grabExpr(
  draw(built$ht, annotation_legend_list = built$legend_list, merge_legend = TRUE)
)

pdf('PCA_and_dotplot.pdf', width = 10, height = 3)
cowplot::plot_grid(PCA_PA, heatmap_grob, labels = c('A','B'))
dev.off()

pdf('PCA_and_dotplot_and_purification.pdf', width = 10, height = 3.5*2)
cowplot::plot_grid(
             cowplot::plot_grid(PCA_PA, heatmap_grob, labels = c('A','B')),
             p_purification_matched, labels = c('','C'), nrow= 2)
dev.off()















PA_HGAP_samples <- embeddings_all_selected %>%
    filter(Cohort %in% c('MINT','Lucas','Sturm NF1'),
           family %in% c('HGAP','PA')) %>%
    pull(Sample_ID) %>% unique()

embeddings_PA_HGAP <- Create_embedding_df(mvalues_combat[,PA_HGAP_samples] , perplexity = 10) %>%
    mutate(Cohort = ifelse(Cohort == 'Sturm','Sturm NF1', Cohort)) %>%
    arrange(Cohort)
pca_full_PA_HGAP<- Run_PCA(mvalues_combat[,PA_HGAP_samples])
var_explained <- (pca_full_PA_HGAP$sdev^2) / sum(pca_full_PA_HGAP$sdev^2) * 100

pdf('PCA_PA_HGAP_HGAP_gliomas_w_batch_correction_family.pdf', height = 3, width = 6)
PCA_PA_HGAP <- embeddings_PA_HGAP %>%
    ggplot(aes(PC1, PC2, fill = family, shape = Cohort, color = Cohort)) +
    geom_point( size = 3,stroke = 1) +
    scale_fill_manual(values = subtype_pal) +
    scale_color_manual(values = cohort_pal) +
    scale_shape_manual(values = c(22,24,21,23)) +
    theme_classic(base_size = 13) +
    labs(
        fill = 'DKFZ class',
        x = paste0('PC1 (',round(var_explained[1]),'%)' ),
        y = paste0('PC2 (',round(var_explained[2]),'%)' ))+
    guides(
        fill = guide_legend(override.aes = list(shape = 21))) +
    geom_point(data = embeddings_PA_HGAP %>% filter(Methylation_Class == 'PA, NF1-associated'),aes(PC1,PC2), shape = 8, color = 'black', inherit.aes = F, size = 1)
PCA_PA_HGAP
dev.off()



pdf('Figure2_tSNE_unedited.pdf', height = 9*0.625 , width = 8*2*0.8)
cowplot::plot_grid(tsne_family, tsne_family_zoom, labels= c('A','B') , nrow = 1, rel_widths = c(1,0.75))
dev.off()





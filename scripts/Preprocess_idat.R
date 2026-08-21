#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Preprocess_idat.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Preprocess .idat files using sesame + mLiftOver.
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# Usage: 
#
# History:
#  11-03-2026: File creation (minfi version)
#  17-08-2026: Rewritten to use sesame + mLiftOver, detection p-value QC,
#              and per-sample platform detection
#              CNV segmentation (sesame::cnSegmentation) added
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(dplyr))
suppressMessages(library(sesame))
suppressMessages(library(sesameData))
suppressMessages(library(BiocParallel))
suppressMessages(library(anndata))

# Use correct reticulate environment
reticulate::use_condaenv(Sys.getenv("CONDA_PREFIX"), required = TRUE)

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if (exists("snakemake")) {
    input <- snakemake@input[[1]]
    Zhou_input <- snakemake@params[['Zhou_probes']]
    CrossReactive_input <- snakemake@params[['CrossReactive_probes']]
    Problematic_input <- snakemake@params[['Problematic_probes']]
    output_adata <- snakemake@output[['adata']]
    threads             <- snakemake@threads
} else {
    input <- '/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/Methylation/samplesheets/Samplesheet_Methylation.csv'
    Zhou_input <- '/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/sslowgrade/workflows/MethylationArray-snake/blocklist/AppendixD_Zhou_et_al_MASKgeneral_list.txt'
    CrossReactive_input <- '/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/sslowgrade/workflows/MethylationArray-snake/blocklist/AppendixE_CrossReactiveProbes_EPICv1.txt'
    Problematic_input <- '/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/sslowgrade/workflows/MethylationArray-snake/blocklist/AppendixF_ProblematicProbes_EPICv1-b5.txt'
    output_adata <- '/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/Methylation//methylation_data.h5ad'
    threads <- 20
}

ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = threads))
register(MulticoreParam(workers = ncores, stop.on.error = FALSE), default = TRUE)

## ---- Prevent CPU oversubscription: MulticoreParam forks `ncores` worker
## processes, but if BLAS/OpenMP inside each worker ALSO spawns its own
## thread pool, you get ncores workers x N threads competing for ncores
## cores - this is what shows up as "way more processes than ncores running"
## in htop/top, and it silently tanks parallel performance. RhpcBLASctl sets
## this at runtime (reliable regardless of BLAS init order); if it isn't
## installed, the env vars below are a decent fallback but only take effect
## if set before the BLAS library initializes.
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(1)
    RhpcBLASctl::omp_set_num_threads(1)
} else {
    Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1", MKL_NUM_THREADS = "1")
}
cat(sprintf(
    "[%s] Using %d worker(s) for parallel steps (parallel::detectCores() reports %d available on this machine)\n",
    Sys.time(), ncores, parallel::detectCores()
))
## ---- Tunable parameters ---------------------------------------------------
DETECTION_PVAL_THRESH  <- 0.05   # pOOBAH detection p-value cutoff. sesame's
                                  # own documented default for pOOBAH is 0.05,
                                  # NOT minfi's usual 0.01 - using 0.01 here
                                  # inflated failure rates substantially.
SAMPLE_OUTLIER_MOD_ZSCORE_THRESH <- 3.5  # a sample is flagged (not dropped)
                                  # if its native detection-failure fraction
                                  # is a modified z-score outlier RELATIVE TO
                                  # ITS OWN PLATFORM's distribution (Iglewicz
                                  # & Hoaglin's standard cutoff). Per-platform
                                  # because EPICv2 has a meaningfully higher
                                  # baseline failure rate than HM450/EPIC even
                                  # in good samples - one global cutoff would
                                  # unfairly flag most of a platform.
PLATFORM_AGE_ORDER     <- c("HM450", "EPIC", "EPICv2")  # oldest to newest -
                                  # TARGET_PLATFORM is deduced from whichever
                                  # of these is OLDEST among the platforms
                                  # actually detected in the samplesheet (2.1b),
                                  # not hardcoded - liftover always goes toward
                                  # the older array, e.g. EPIC+EPICv2 -> EPIC.
DROP_SEX_CHROM         <- TRUE
PROBE_BLACKLIST_FAIL_FRACTION <- 0.1   # beyond the curated blacklists, also
                                        # permanently exclude any probe that
                                        # FAILED DETECTION in more than this
                                        # fraction of RELIABLE samples (i.e.
                                        # excluding samples flagged as
                                        # detection-failure outliers - see
                                        # SAMPLE_OUTLIER_MOD_ZSCORE_THRESH).
                                        # Decoupled from any missingness
                                        # tolerance - a badly-behaved probe
                                        # stays excluded here regardless.
CHUNK_SIZE             <- 100   # per-platform batch size for openSesame
                                  # calls, kept modest for memory safety
RUN_CNV                <- TRUE   # run sesame::cnSegmentation() per sample
CNV_TILEWIDTH           <- 50000 # bin width (bp) for cnSegmentation, sesame's own default
 
#-------------------------------------------------------------------------------
# 1.1 Read sample sheet (only the samples listed here are ever touched)
#-------------------------------------------------------------------------------
samplesheet <- read.delim(input, sep = ',') %>%
    mutate(
        idat_basename = gsub("_Red\\.idat$", "", idat_red),
        Basename_full = idat_basename,
        batch         = basename(dirname(idat_red))
    )
 
stopifnot(!any(duplicated(samplesheet$sample)))
stopifnot(all(file.exists(paste0(samplesheet$Basename_full, "_Grn.idat"))))
cat(sprintf("[%s] Sample sheet loaded: %d samples\n", Sys.time(), nrow(samplesheet)))
 
#-------------------------------------------------------------------------------
# 1.2 Probe blacklists
#-------------------------------------------------------------------------------
# NB: the Zhou et al. MASKgeneral list already folds in SNP-affected probes,
# so there is no separate dropLociWithSnps()-style step here (sesame has no
# direct equivalent, and it would be redundant with this blacklist).
Zhou_probes          <- read.delim(Zhou_input, col.names = 'Probe', header = FALSE)
CrossReactive_probes <- read.delim(CrossReactive_input)
Problematic_probes   <- read.delim(Problematic_input, col.names = 'Probe')
Filter_probes <- unique(c(Zhou_probes$Probe, CrossReactive_probes$Probe, Problematic_probes$Probe))
cat(sprintf("[%s] Blacklist probes (Zhou + cross-reactive + problematic): %d\n",
            Sys.time(), length(Filter_probes)))
 
#-------------------------------------------------------------------------------
# 2.1 Per-sample platform detection
#-------------------------------------------------------------------------------
detect_platform <- function(prefix) {
    sdf <- tryCatch(sesame::readIDATpair(prefix), error = function(e) NULL)
    if (is.null(sdf)) return(NA_character_)
    plat <- attr(sdf, "platform")
    rm(sdf); gc(FALSE)
    if (is.null(plat) || length(plat) == 0) return(NA_character_)
    as.character(plat)
}
 
samplesheet$Platform <- unlist(bplapply(samplesheet$Basename_full, detect_platform, BPPARAM = bpparam()))
samplesheet <- samplesheet %>% filter(Platform %in% c("HM450", "EPIC", "EPICv2"))
 
#-------------------------------------------------------------------------------
# 2.1b Deduce TARGET_PLATFORM from what's actually in the samplesheet, then
#      build the target manifest. Liftover always goes toward the OLDEST
#      platform present (per PLATFORM_AGE_ORDER) - e.g. HM450+EPICv2 -> HM450,
#      EPIC+EPICv2 (no HM450) -> EPIC, EPICv2-only -> EPICv2 (no liftover
#      needed at all in that case). This has to happen AFTER platform
#      detection (2.1), since it depends on what was actually found.
#-------------------------------------------------------------------------------
detected_platforms <- unique(samplesheet$Platform)
TARGET_PLATFORM <- PLATFORM_AGE_ORDER[PLATFORM_AGE_ORDER %in% detected_platforms][1]
target_anno       <- as.data.frame(sesameData_getManifestGRanges(TARGET_PLATFORM))
target_probe_ids  <- rownames(target_anno)
sexchrom_probes   <- rownames(target_anno)[target_anno$seqnames %in% c("chrX", "chrY")]
cat(sprintf("[%s] %s sex-chromosome probes: %d (DROP_SEX_CHROM=%s)\n",
            Sys.time(), TARGET_PLATFORM, length(sexchrom_probes), DROP_SEX_CHROM))
 
 
#-------------------------------------------------------------------------------
# 2.2 Per-platform sesame processing: plain QCDPB betas, lifted onto the
#     common TARGET_PLATFORM space (2.1b), plus a lifted detection-failure
#     indicator matrix. Detection p-values (pOOBAH) drive two independent
#     things:
#       - sample-level QC in 2.3 (native, pre-liftover fail fraction)
#       - the blacklist extension in 2.4 (post-liftover, per-probe fail
#         fraction across kept samples)
#     They never touch the beta values themselves - betas stay unmasked.
#-------------------------------------------------------------------------------
process_chunk <- function(prefixes, sample_ids, platform) {
 
    ## ---- Betas: full QCDPB pipeline, unmasked
    betas_native <- openSesame(
        prefixes, prep = "QCDPB",
        func = function(sdf) sesame::getBetas(sdf, mask = FALSE),
        BPPARAM = bpparam()
    )
    if (is.null(dim(betas_native))) {
        betas_native <- matrix(betas_native, ncol = 1, dimnames = list(names(betas_native), sample_ids))
    } else {
        colnames(betas_native) <- sample_ids
    }
 
    ## ---- Detection p-values: QCD only (no P-masking baked in), pOOBAH called
    ## explicitly so we get the raw per-probe p-values rather than a mask.
    pvals_native <- openSesame(
        prefixes, prep = "QCD",
        func = function(sdf) sesame::pOOBAH(sdf, return.pval = TRUE),
        BPPARAM = bpparam()
    )
    if (is.null(dim(pvals_native))) {
        pvals_native <- matrix(pvals_native, ncol = 1, dimnames = list(names(pvals_native), sample_ids))
    } else {
        colnames(pvals_native) <- sample_ids
    }
 
    # A probe FAILS detection when p > threshold; NA also counts as failed
    # since detection genuinely couldn't be assessed for it.
    failed_native <- (pvals_native > DETECTION_PVAL_THRESH) | is.na(pvals_native)
    failed_native <- matrix(as.numeric(failed_native), nrow = nrow(failed_native),
                             dimnames = dimnames(failed_native))
 
    n_probes_native <- nrow(betas_native)
    n_failed_native  <- colSums(failed_native, na.rm = TRUE)
 
    ## ---- Lift betas AND the failure indicator onto the common target-
    ## platform space. Same row-alignment logic for both: mLiftOver doesn't
    ## guarantee a row for every target probe, so explicitly reindex its output.
    if (identical(platform, TARGET_PLATFORM)) {
        betas_target  <- betas_native[match(target_probe_ids, rownames(betas_native)), , drop = FALSE]
        failed_target <- failed_native[match(target_probe_ids, rownames(failed_native)), , drop = FALSE]
    } else {
        betas_target_raw  <- mLiftOver(betas_native, TARGET_PLATFORM, impute = FALSE, BPPARAM = bpparam())
        failed_target_raw <- mLiftOver(failed_native, TARGET_PLATFORM, impute = FALSE, BPPARAM = bpparam())
        betas_target  <- betas_target_raw[match(target_probe_ids, rownames(betas_target_raw)), , drop = FALSE]
        failed_target <- failed_target_raw[match(target_probe_ids, rownames(failed_target_raw)), , drop = FALSE]
        rm(betas_target_raw, failed_target_raw)
    }
    rownames(betas_target)  <- target_probe_ids
    rownames(failed_target) <- target_probe_ids
    stopifnot(identical(rownames(betas_target), target_probe_ids))
 
    list(
        betas_native  = betas_native,
        betas_target  = betas_target,
        failed_target = failed_target,
        qc = data.frame(
            sample = sample_ids, Platform = platform,
            n_probes_native = n_probes_native,
            n_failed_native  = as.integer(n_failed_native),
            frac_failed_native = as.numeric(n_failed_native) / n_probes_native,
            row.names = NULL
        )
    )
}
 
beta_native_chunks <- list()   # keyed by platform - kept in each platform's OWN native probe space
beta_target_chunks   <- list()
failed_target_chunks <- list()
sample_qc <- data.frame(
    sample = character(), Platform = character(),
    n_probes_native = integer(), n_failed_native = integer(),
    frac_failed_native = double(), stringsAsFactors = FALSE
)

for (plat in c("HM450", "EPIC", "EPICv2")) {
    plat_sheet <- samplesheet %>% filter(Platform == plat) 
    if (nrow(plat_sheet) == 0) next
    n_chunks <- ceiling(nrow(plat_sheet) / CHUNK_SIZE)
    cat(sprintf("[%s] Platform %s: %d samples in %d chunk(s) of <=%d\n",
                Sys.time(), plat, nrow(plat_sheet), n_chunks, CHUNK_SIZE))
    plat_native_chunks <- list()
    for (i in seq_len(n_chunks)) {
        idx <- ((i - 1) * CHUNK_SIZE + 1):min(i * CHUNK_SIZE, nrow(plat_sheet))
        res <- process_chunk(
            prefixes   = plat_sheet$Basename_full[idx],
            sample_ids = plat_sheet$sample[idx],
            platform   = plat
        )
        key <- paste0(plat, "_", i)
        plat_native_chunks[[key]]    <- res$betas_native
        beta_target_chunks[[key]]    <- res$betas_target
        failed_target_chunks[[key]]  <- res$failed_target
        sample_qc <- rbind(sample_qc, res$qc)
        rm(res); gc(FALSE)
    }
    # native probe IDs are identical across chunks of the same platform
    # (same manifest), so a plain cbind is safe here.
    beta_native_chunks[[plat]] <- do.call(cbind, plat_native_chunks)
    rm(plat_native_chunks); gc(FALSE)
}
 
beta_target_all   <- do.call(cbind, beta_target_chunks)
failed_target_all <- do.call(cbind, failed_target_chunks)
# restore original samplesheet row order
beta_target_all   <- beta_target_all[,   samplesheet$sample, drop = FALSE]
failed_target_all <- failed_target_all[, samplesheet$sample, drop = FALSE]
rm(beta_target_chunks, failed_target_chunks); gc(FALSE)
 
#-------------------------------------------------------------------------------
# 2.3 Sample-level QC: FLAG (not drop) samples whose native detection-
#     failure fraction is a strong outlier relative to their OWN platform's
#     distribution. No sample is removed here - every sample that had a
#     detectable platform stays through to the final AnnData object; the
#     flag and the numbers behind it are attached to `obs` instead, same
#     "annotate rather than subset" approach used for probes in 2.5/2.6.
#-------------------------------------------------------------------------------
modified_zscore <- function(x) {
    med   <- median(x, na.rm = TRUE)
    mad_x <- mad(x, center = med, constant = 1, na.rm = TRUE)  # raw MAD (no normal-consistency scaling)
    if (mad_x == 0) return(rep(0, length(x)))  # no spread within this platform to compare against
    0.6745 * (x - med) / mad_x
}
 
sample_qc <- sample_qc %>%
    group_by(Platform) %>%
    mutate(
        platform_n                   = n(),
        platform_median_frac_failed  = median(frac_failed_native),
        modified_zscore              = modified_zscore(frac_failed_native),
        # Only flag samples that fail MORE than their platform's typical
        # rate - an unusually clean sample isn't a QC problem.
        sample_qc_outlier            = modified_zscore > SAMPLE_OUTLIER_MOD_ZSCORE_THRESH
    ) %>%
    ungroup() %>%
    as.data.frame()
 
 
cat(sprintf(
    "[%s] %d/%d samples flagged as detection-failure outliers (modified z-score > %.1f within their own platform)\n",
    Sys.time(), sum(sample_qc$sample_qc_outlier), nrow(sample_qc), SAMPLE_OUTLIER_MOD_ZSCORE_THRESH
))
if (any(sample_qc$sample_qc_outlier)) {
    print(sample_qc %>% filter(sample_qc_outlier) %>%
              select(sample, Platform, frac_failed_native, platform_median_frac_failed, modified_zscore))
}
 
all_samples <- samplesheet$sample  # every sample stays - nothing dropped here
 
#-------------------------------------------------------------------------------
# 2.4 Extend the probe blacklist based on DETECTION-P-VALUE failure rate
#-------------------------------------------------------------------------------
# "Failed" = pOOBAH detection p-value > DETECTION_PVAL_THRESH, NOT the
# generic NA pattern from mLiftOver structural gaps (that says nothing
# about detection quality). Beyond the curated Zhou/cross-reactive/
# problematic blacklist from 1.2, permanently exclude any probe that failed
# DETECTION in more than PROBE_BLACKLIST_FAIL_FRACTION of ALL samples -
# sample_qc_outlier (2.3) is purely a downstream annotation here and isn't
# used to filter this denominator, to keep the two QC layers independent
# and simple. This is deliberately independent of any missingness tolerance
# you might apply downstream - a badly-behaved probe stays excluded here
# regardless. Reassigns `Filter_probes` in place so every later use of it
# picks up the extension automatically.
frac_failed_detection_per_probe <- rowSums(failed_target_all) / ncol(failed_target_all)
 
newly_blacklisted <- names(which(frac_failed_detection_per_probe > PROBE_BLACKLIST_FAIL_FRACTION))
n_new_blacklisted <- length(setdiff(newly_blacklisted, Filter_probes))
Filter_probes <- union(Filter_probes, newly_blacklisted)
 
cat(sprintf(
    "[%s] Extended blacklist: %d probes newly added (failed detection in >%.1f%% of %d samples). Total blacklist now: %d probes.\n",
    Sys.time(), n_new_blacklisted, 100 * PROBE_BLACKLIST_FAIL_FRACTION,
    length(all_samples), length(Filter_probes)
))
 
#-------------------------------------------------------------------------------
# 2.5 Probe-level QC flags - computed for EVERY probe in TARGET_PLATFORM's
#     manifest, nothing is subset out of the matrix. `passed_qc` is the
#     combined recommendation (curated + detection-failure blacklist,
#     sex-chrom, never-measured); everything it's built from is also kept
#     as its own column so you can recombine the criteria differently later
#     without re-running the script.
#-------------------------------------------------------------------------------
frac_missing_per_probe <- rowMeans(is.na(beta_target_all))
never_measured_probes  <- names(which(frac_missing_per_probe >= 1))
if (length(never_measured_probes) > 0) {
    cat(sprintf(
        "[%s] %d probe(s) have no measurement in any kept sample (structural liftover gap) - flagged, not removed\n",
        Sys.time(), length(never_measured_probes)
    ))
}
 
curated_blacklist <- unique(c(Zhou_probes$Probe, CrossReactive_probes$Probe, Problematic_probes$Probe))
 
#-------------------------------------------------------------------------------
# 2.6 Build var (probe metadata): manifest annotation + QC flag columns,
#     for the full TARGET_PLATFORM probe set - every probe stays in the object.
#-------------------------------------------------------------------------------
probe_metadata <- target_anno %>%
    mutate(across(everything(), ~ ifelse(is.na(.), "", as.character(.))))
 
probe_metadata$frac_missing          <- frac_missing_per_probe[rownames(probe_metadata)]
probe_metadata$never_measured        <- rownames(probe_metadata) %in% never_measured_probes
probe_metadata$frac_failed_detection <- frac_failed_detection_per_probe[rownames(probe_metadata)]
probe_metadata$on_curated_blacklist  <- rownames(probe_metadata) %in% curated_blacklist
probe_metadata$on_extended_blacklist <- rownames(probe_metadata) %in% Filter_probes  # curated + detection-failure additions
probe_metadata$sex_chrom             <- rownames(probe_metadata) %in% sexchrom_probes
probe_metadata$passed_qc <- !probe_metadata$on_extended_blacklist &
                             !probe_metadata$never_measured &
                             !(DROP_SEX_CHROM & probe_metadata$sex_chrom)
 
cat(sprintf(
    "[%s] %d/%d %s probes pass QC (flagged in var$passed_qc, none removed from the object)\n",
    Sys.time(), sum(probe_metadata$passed_qc), nrow(probe_metadata), TARGET_PLATFORM
))
 
#-------------------------------------------------------------------------------
# 2.7 Assemble beta matrices for storage - beta values only throughout
#     (M-values are a one-line transform of beta if/when you need them:
#     log2(beta / (1 - beta)), so there's no need to also store them).
#-------------------------------------------------------------------------------
# X: the fully normalized, harmonized beta matrix for every TARGET_PLATFORM
# probe, BEFORE any probe-level filtering ("pre-blacklisted") - this is the
# always-recoverable "normalized data with all probes" view.
beta_prelift_blacklist <- beta_target_all
 
#-------------------------------------------------------------------------------
# 2.8 CNV segmentation (sesame::cnSegmentation), per sample, against a
#     platform-appropriate normal reference. Runs off a FRESH read of each
#     sample's IDATs (raw SigDF, no QCDPB prep) rather than reusing the
#     already-normalized signal from 2.2 - cnSegmentation does its own
#     normalization against the normal reference set, and pre-processing
#     (noob background correction, dye-bias correction) would distort the
#     total-intensity signal that copy-number calling depends on. This is
#     therefore a genuine extra I/O pass over every IDAT pair.
#
# Normal reference: sesame ships built-in defaults for EPIC and EPICv2
# (cnv_normal_default() below, taken from sesame's own source) but NOT for
# HM450 - sesameData does publish TCGA-derived HM450 normal-tissue signal
# sets, so those are tried as a fallback. Bioconductor/sesameData dataset
# names have shifted across sesame versions (SigSet -> SigDF), so this
# tries a few known candidate names and uses whichever your installed
# sesameData actually has; if none match, it errors with a pointer to
# `sesameData::sesameDataList()` so you can find the current name and set
# it explicitly.
#
# There's no Rds output in this script anymore, so the raw CNSegment
# objects (the only thing sesame::visualizeSegments() accepts) are NOT
# preserved anywhere - only flattened data.frame/matrix versions make it
# into adata$uns, same as everything else here.
#-------------------------------------------------------------------------------
if (RUN_CNV) {
 
cnv_normal_default <- function(platform) {
    candidates <- switch(platform,
        "EPICv2" = list(list(name = "EPICv2.8.SigDF",
                              select = c("GM12878_206909630042_R08C01", "GM12878_206909630040_R03C01"))),
        "EPIC"   = list(list(name = "EPIC.5.SigDF.normal", select = NULL),
                        list(name = "EPIC.5.normal",        select = NULL)),
        "HM450"  = list(list(name = "HM450.10.SigDF.normal",     select = NULL),
                        list(name = "HM450.10.TCGA.PAAD.normal", select = NULL),
                        list(name = "HM450.10.TCGA.BLCA.normal", select = NULL)),
        stop(sprintf("No known built-in CNV normal reference for platform %s", platform))
    )
    for (cand in candidates) {
        sdfs.normal <- tryCatch(sesameDataGet(cand$name), error = function(e) NULL)
        if (is.null(sdfs.normal)) next
        if (!is.null(cand$select)) sdfs.normal <- sdfs.normal[cand$select]
        # Older/alternately-shaped sesameData objects sometimes wrap the
        # signal sets under a $sset or $ssets field rather than being the
        # list of SigDF/SigSet objects directly - handle both shapes.
        if (is.list(sdfs.normal) && !is.null(sdfs.normal$ssets)) return(sdfs.normal$ssets)
        if (is.list(sdfs.normal) && !is.null(sdfs.normal$sset))  return(list(sdfs.normal$sset))
        return(sdfs.normal)
    }
    stop(sprintf(
        paste0("Could not find a built-in CNV normal reference for platform %s in your installed ",
               "sesameData. Run sesameData::sesameDataList() to see what's available and update ",
               "cnv_normal_default() with the current dataset name."),
        platform
    ))
}

run_cnv_for_sample <- function(prefix, sample_id, sdfs.normal) {
    tryCatch({
        sdf <- sesame::readIDATpair(prefix)        
        seg <- sesame::cnSegmentation(sdf, sdfs.normal = sdfs.normal, tilewidth = CNV_TILEWIDTH)
        list(sample = sample_id, seg = seg, error = NA_character_)
    }, error = function(e) {
        list(sample = sample_id, seg = NULL, error = conditionMessage(e))
    })
}

 
cnv_results_by_platform <- list()
for (plat in unique(samplesheet$Platform)) {
    plat_sheet <- samplesheet %>% filter(Platform == plat)
    cat(sprintf("[%s] CNV segmentation: %s, %d sample(s)\n", Sys.time(), plat, nrow(plat_sheet)))
 
    sdfs.normal <- tryCatch(cnv_normal_default(plat), error = function(e) {
        warning(sprintf("Skipping CNV for platform %s: %s", plat, conditionMessage(e)))
        NULL
    })
    if (is.null(sdfs.normal)) next
 
    cnv_results_by_platform[[plat]] <- bplapply(
        seq_len(nrow(plat_sheet)),
        function(i) run_cnv_for_sample(plat_sheet$Basename_full[i], plat_sheet$sample[i], sdfs.normal),
        BPPARAM = bpparam()
    )
}
    
    cnv_results <- unlist(cnv_results_by_platform, recursive = FALSE)

    
cnv_failures <- Filter(function(r) !is.na(r$error), cnv_results)
if (length(cnv_failures) > 0) {
    warning(sprintf("CNV segmentation failed for %d/%d sample(s): %s",
                     length(cnv_failures), length(cnv_results),
                     paste(sapply(cnv_failures, `[[`, "sample"), collapse = ", ")))
}
cnv_ok <- Filter(function(r) is.na(r$error), cnv_results)
cat(sprintf("[%s] CNV segmentation succeeded for %d/%d samples\n",
            Sys.time(), length(cnv_ok), length(cnv_results)))
 
# ---- Flatten into h5ad-friendly structures (plain data.frames/matrices).
# NB: field names below (bin.coords/bin.signals/seg.signals) match sesame's
# documented CNSegment structure as of this writing - if your installed
# sesame version differs, `str(cnv_ok[[1]]$seg)` will show the actual
# fields to adjust this against.
cnv_bin_signals_by_platform <- list()
cnv_bin_coords_by_platform  <- list()
cnv_segments_list <- list()
for (r in cnv_ok) {
    seg <- r$seg
    plat <- samplesheet$Platform[match(r$sample, samplesheet$sample)]
 
    bin_signal <- tryCatch(as.numeric(seg$bin.signals), error = function(e) NULL)
    if (!is.null(bin_signal)) {
        if (is.null(cnv_bin_signals_by_platform[[plat]])) cnv_bin_signals_by_platform[[plat]] <- list()
        cnv_bin_signals_by_platform[[plat]][[r$sample]] <- bin_signal
    }
    if (is.null(cnv_bin_coords_by_platform[[plat]])) {
        cnv_bin_coords_by_platform[[plat]] <- tryCatch(as.data.frame(seg$bin.coords), error = function(e) NULL)
    }
 
    seg_df <- tryCatch(as.data.frame(seg$seg.signals), error = function(e) NULL)
    if (!is.null(seg_df)) {
        seg_df$sample   <- r$sample
        seg_df$Platform <- plat
        cnv_segments_list[[r$sample]] <- seg_df
    }
}
 
cnv_bin_signals_by_platform <- lapply(cnv_bin_signals_by_platform, function(sample_list) {
    m <- do.call(cbind, sample_list)
    colnames(m) <- names(sample_list)
    as.data.frame(t(m)) %>% tibble::rownames_to_column("sample")
})
cnv_segments <- if (length(cnv_segments_list) > 0) do.call(rbind, cnv_segments_list) else NULL
} # if (RUN_CNV)
#-------------------------------------------------------------------------------
# 3.0 Create AnnData object
#-------------------------------------------------------------------------------
# Attach the sample-level QC numbers/flag from 2.3 to obs - every sample
# stays (this is a left_join, not a filter), so obs carries the same rows
# as before plus n_probes_native, n_failed_native, frac_failed_native,
# platform_median_frac_failed, modified_zscore, and sample_qc_outlier.

samplesheet <- samplesheet %>%
    left_join(
        sample_qc %>% select(sample, n_probes_native, n_failed_native, frac_failed_native,
                              platform_median_frac_failed, modified_zscore, sample_qc_outlier),
        by = "sample"
    )

adata <- anndata::AnnData(
    X   = t(beta_prelift_blacklist),   # normalized beta, all TARGET_PLATFORM probes, pre-blacklist
    obs = samplesheet,
    var = probe_metadata,
    )

rownames(adata$obs) <- samplesheet$sample


# Native, pre-liftover beta values - kept in `uns` (not `layers`) because
# each platform has its own probe space, which can't share var with X:
# AnnData requires every layer to match X's shape exactly, and a 450K/EPIC/
# EPICv2 mix genuinely doesn't share one probe axis pre-liftover. One
# data.frame per platform, each with its own native probe IDs as a column.
adata$uns[['beta_prelift_by_platform']] <- lapply(beta_native_chunks, function(m) {
    data.frame(m, check.names = FALSE) %>% tibble::rownames_to_column("Probe_ID")
})
 
if (RUN_CNV) {
    # bin.coords/bin.signals differ by platform (different probe coverage
    # feeds each genomic bin), so these stay split per platform, same
    # reasoning as beta_prelift_by_platform above. cnv_segments is one long
    # (tidy) data.frame across all samples/platforms, since segment COUNT
    # varies per sample and can't be a fixed-shape matrix.
    adata$uns[['cnv_bin_coords_by_platform']]  <- cnv_bin_coords_by_platform
    adata$uns[['cnv_bin_signals_by_platform']] <- cnv_bin_signals_by_platform
    adata$uns[['cnv_segments']] <- cnv_segments
    if (length(cnv_failures) > 0) {
        adata$uns[['cnv_failures']] <- data.frame(
            sample = sapply(cnv_failures, `[[`, "sample"),
            error  = sapply(cnv_failures, `[[`, "error")
        )
    }
}


#-------------------------------------------------------------------------------
# 4.0 Write outputs
#-------------------------------------------------------------------------------
write_h5ad(adata, output_adata)

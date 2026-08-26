#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Preprocess_idat.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Preprocess .idat files listed in a samplesheet using sesame + mLiftOver.
#
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# Usage: invoked by the Snakemake rule `Preprocess_idat`
#
# History:
#  11-03-2026: File creation (minfi version)
#  17-08-2026: Rewritten to use sesame + mLiftOver, detection p-value QC,
#              and per-sample platform detection
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

library(minfi)
purified <- read.metharray('~/mnt/BIGR_home/MINT/data/idat/purified/purified/MINT_M06_pur05/206467010175_R05C01')
nonpurified <- read.metharray('~/mnt/BIGR_home/MINT/data/idat/purified/unclass_unpurified/MINT_M06/206467010175_R05C01')
lost <- setdiff(rownames(nonpurified), rownames(purified))
length(lost)   # expect 404

ctrl <- minfi::getProbeInfo(nonpurified, type = "Control")
sum(lost %in% ctrl$Address)          # any control probes lost?
table(ctrl$Type[ctrl$Address %in% lost])


mset <- preprocessIllumina(purified, bg.correct = TRUE, normalize = "controls")
b    <- getBeta(mset, offset = 100)
summary(b)
mean(is.na(b))
mset_o <- preprocessIllumina(nonpurified, bg.correct = TRUE, normalize = "controls")
b_o <- getBeta(mset_o, offset = 100)
shared <- intersect(rownames(b_o), rownames(b))

plot(density(b_o[shared, 1]), main = "beta", lwd = 2)
lines(density(b[shared, 1]), col = "red", lwd = 2)
legend("topright", c("original", "purified"), col = c("black", "red"), lwd = 2)

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if (exists("snakemake")) {
    input <- snakemake@input[[1]]
    Zhou_input <- snakemake@params[['Zhou_probes']]
    CrossReactive_input <- snakemake@params[['CrossReactive_probes']]
    Problematic_input <- snakemake@params[['Problematic_probes']]
    output_adata <- snakemake@output[['adata']]
    output_Mset <- snakemake@output[['Mset']]
    threads <- snakemake@threads
} else {
    input               <- '/home/jurriaan/Projects/MINT/data/samplesheets/samplesheet_methylation.csv'
    Zhou_input          <- '/data/Resources/EPIC/manifest/AppendixD_Zhou_et_al_MASKgeneral_list.txt'
    CrossReactive_input <- '/data/Resources/EPIC/manifest/AppendixE_CrossReactiveProbes_EPICv1.txt'
    Problematic_input   <- '/data/Resources/EPIC/manifest/AppendixF_ProblematicProbes_EPICv1-b5.txt'
    output_adata        <- 'output/methylation/methylation_data.h5ad'
    output_Mset         <- 'output/methylation/methylation_object.Rds'
    threads             <- 2
}

ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = threads))
register(MulticoreParam(workers = ncores, stop.on.error = FALSE), default = TRUE)

## ---- Tunable parameters ---------------------------------------------------
DETECTION_PVAL_THRESH  <- 0.01   # pOOBAH detection p-value cutoff
MAX_SAMPLE_MISSING     <- 0.10   # exclude a sample if >10% of its NATIVE
                                  # probes fail detection (same threshold
                                  # the old minfi script actually used,
                                  # despite its comment saying 5%)
MAX_PROBE_FAIL_FRACTION <- 0     # a probe is kept only if it fails
                                  # detection in at most this fraction of
                                  # KEPT samples; 0 = must pass detection
                                  # in every kept sample (matches the old
                                  # script's `== ncol(...)` behaviour)
TARGET_PLATFORM        <- "HM450"
DROP_SEX_CHROM         <- TRUE
CHUNK_SIZE             <- 500    # per-platform batch size for openSesame
                                  # calls, kept modest for memory safety

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
# 1.2 Probe blacklists + reference (HM450) manifest
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

hm450_anno       <- as.data.frame(sesameData_getManifestGRanges(TARGET_PLATFORM))
hm450_probe_ids  <- rownames(hm450_anno)
sexchrom_probes  <- rownames(hm450_anno)[hm450_anno$seqnames %in% c("chrX", "chrY")]
cat(sprintf("[%s] HM450 sex-chromosome probes: %d (DROP_SEX_CHROM=%s)\n",
            Sys.time(), length(sexchrom_probes), DROP_SEX_CHROM))

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

cat(sprintf("[%s] Detecting array platform for %d samples...\n", Sys.time(), nrow(samplesheet)))
samplesheet$Platform <- unlist(bplapply(samplesheet$Basename_full, detect_platform, BPPARAM = bpparam()))

n_undetected <- sum(is.na(samplesheet$Platform))
if (n_undetected > 0) {
    warning(sprintf(
        "%d/%d samples had no platform attribute from sesame - falling back to probe-count classification for those only.",
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

cat("[Platform distribution]\n")
print(table(samplesheet$Platform, useNA = "ifany"))

n_before <- nrow(samplesheet)
samplesheet <- samplesheet %>% filter(Platform %in% c("HM450", "EPIC", "EPICv2"))
if (nrow(samplesheet) < n_before) {
    warning(sprintf("Dropped %d sample(s) with an unrecognized/undetectable platform.",
                     n_before - nrow(samplesheet)))
}
stopifnot(nrow(samplesheet) > 0)

#-------------------------------------------------------------------------------
# 2.2 Per-platform sesame processing: betas (QCDPB, mask=FALSE) + detection
#     p-values (QCD + pOOBAH), both lifted onto the common HM450 space
#-------------------------------------------------------------------------------
process_chunk <- function(prefixes, sample_ids, platform) {

    ## ---- Betas: full QCDPB pipeline, unmasked so liftover has data to work with
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
    n_failed_native  <- colSums(failed_native)

    ## ---- Lift both matrices onto the common HM450 probe space
    if (identical(platform, TARGET_PLATFORM)) {
        betas_hm450  <- betas_native[match(hm450_probe_ids, rownames(betas_native)), , drop = FALSE]
        failed_hm450 <- failed_native[match(hm450_probe_ids, rownames(failed_native)), , drop = FALSE]
    } else {
        betas_hm450_raw  <- mLiftOver(betas_native, TARGET_PLATFORM, impute = FALSE, BPPARAM = bpparam())
        failed_hm450_raw <- mLiftOver(failed_native, TARGET_PLATFORM, impute = FALSE, BPPARAM = bpparam())
        betas_hm450  <- betas_hm450_raw[match(hm450_probe_ids, rownames(betas_hm450_raw)), , drop = FALSE]
        failed_hm450 <- failed_hm450_raw[match(hm450_probe_ids, rownames(failed_hm450_raw)), , drop = FALSE]
        rm(betas_hm450_raw, failed_hm450_raw)
    }
    rownames(betas_hm450)  <- hm450_probe_ids
    rownames(failed_hm450) <- hm450_probe_ids
    stopifnot(identical(rownames(betas_hm450), hm450_probe_ids))

    list(
        betas_hm450  = betas_hm450,
        failed_hm450 = failed_hm450,
        qc = data.frame(
            sample = sample_ids, Platform = platform,
            n_probes_native = n_probes_native,
            n_failed_native  = as.integer(n_failed_native),
            frac_failed_native = as.numeric(n_failed_native) / n_probes_native,
            row.names = NULL
        )
    )
}

beta_hm450_chunks   <- list()
failed_hm450_chunks <- list()
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

    for (i in seq_len(n_chunks)) {
        idx <- ((i - 1) * CHUNK_SIZE + 1):min(i * CHUNK_SIZE, nrow(plat_sheet))
        res <- process_chunk(
            prefixes   = plat_sheet$Basename_full[idx],
            sample_ids = plat_sheet$sample[idx],
            platform   = plat
        )
        key <- paste0(plat, "_", i)
        beta_hm450_chunks[[key]]   <- res$betas_hm450
        failed_hm450_chunks[[key]] <- res$failed_hm450
        sample_qc <- rbind(sample_qc, res$qc)
        rm(res); gc(FALSE)
    }
}

beta_hm450_all   <- do.call(cbind, beta_hm450_chunks)
failed_hm450_all <- do.call(cbind, failed_hm450_chunks)
# restore original samplesheet row order
beta_hm450_all   <- beta_hm450_all[,   samplesheet$sample, drop = FALSE]
failed_hm450_all <- failed_hm450_all[, samplesheet$sample, drop = FALSE]
rm(beta_hm450_chunks, failed_hm450_chunks); gc(FALSE)

#-------------------------------------------------------------------------------
# 2.3 Sample-level QC: drop samples with too many failed NATIVE probes
#-------------------------------------------------------------------------------
failed_samples <- sample_qc$sample[sample_qc$frac_failed_native > MAX_SAMPLE_MISSING]

if (length(failed_samples) > 0) {
    message("Failed samples: ", paste(failed_samples, collapse = ", "))
}

kept_samples <- setdiff(samplesheet$sample, failed_samples)
stopifnot(length(kept_samples) > 0)

beta_hm450_all   <- beta_hm450_all[,   kept_samples, drop = FALSE]
failed_hm450_all <- failed_hm450_all[, kept_samples, drop = FALSE]
samplesheet <- samplesheet %>% filter(sample %in% kept_samples)
# keep sample QC only for samples we kept, indexed for later reporting
sample_qc <- sample_qc %>% filter(sample %in% kept_samples)

cat(sprintf("[%s] %d/%d samples kept after native detection-failure QC (>%.0f%% cutoff)\n",
            Sys.time(), length(kept_samples), length(kept_samples) + length(failed_samples),
            100 * MAX_SAMPLE_MISSING))

#-------------------------------------------------------------------------------
# 2.4 Probe-level QC: keep probes that pass detection in (nearly) all kept
#     samples, then apply the curated blacklist + sex-chromosome removal
#-------------------------------------------------------------------------------
frac_failed_detection_per_probe <- rowSums(failed_hm450_all) / ncol(failed_hm450_all)

keep_probes <- names(which(frac_failed_detection_per_probe <= MAX_PROBE_FAIL_FRACTION))
keep_probes <- setdiff(keep_probes, Filter_probes)
if (DROP_SEX_CHROM) keep_probes <- setdiff(keep_probes, sexchrom_probes)
keep_probes <- sort(keep_probes)

cat(sprintf("[%s] Remaining probes after filtering: %d\n", Sys.time(), length(keep_probes)))
stopifnot(length(keep_probes) > 0)

#-------------------------------------------------------------------------------
# 2.5 Extract final data
#-------------------------------------------------------------------------------
beta_to_mvalue <- function(beta, offset = 1e-6) {
    beta_clipped <- pmin(pmax(beta, offset), 1 - offset)
    log2(beta_clipped / (1 - beta_clipped))
}

beta_values <- beta_hm450_all[keep_probes, , drop = FALSE]
m_values    <- beta_to_mvalue(beta_values)

# raw_beta_values: the unfiltered, mask=FALSE beta matrix in the common
# HM450 space (i.e. before blacklist/detection/sex-chrom filtering) - the
# harmonized analogue of getBeta(raw_intensity_data) in the old script.
raw_beta_values <- beta_hm450_all

# probe metadata (var), restricted to the final kept probe set
probe_metadata <- hm450_anno[keep_probes, , drop = FALSE] %>%
    mutate(across(everything(), ~ ifelse(is.na(.), "", as.character(.))))

#-------------------------------------------------------------------------------
# 3.0 Create AnnData object
#-------------------------------------------------------------------------------
print(samplesheet)
adata <- anndata::AnnData(
    X   = t(m_values),
    obs = samplesheet,
    var = probe_metadata,
)

adata$uns[['beta_raw']] <- data.frame(raw_beta_values) %>% tibble::rownames_to_column()
adata$layers[['beta']] <- t(beta_values)

#-------------------------------------------------------------------------------
# 4.0 Write outputs
#-------------------------------------------------------------------------------
write_h5ad(adata, output_adata)

# There is no sesame equivalent of minfi's MethylSet class, so the "Mset"
# output is now a plain list carrying everything a downstream script would
# need: final filtered beta/M-value matrices, probe metadata, and the
# post-QC sample sheet (with detected Platform).
methylation_object <- list(
    beta            = beta_values,
    mvalue          = m_values,
    probe_metadata  = probe_metadata,
    samplesheet     = samplesheet,
    sample_qc       = sample_qc
)
saveRDS(methylation_object, output_Mset)

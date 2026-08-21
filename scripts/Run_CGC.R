#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Run_CGC.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Continuous Grading Coefficient (CGCpsi) for astrocytoma, computed directly
# from the preprocessed methylation AnnData rather than from raw idats.
# Adapted from CGC-Psi.R,
# https://github.com/ErasmusMC-Neuro-Oncology/Continuous_Grading_Classifier
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv: envs/CGC.yaml
# Usage: snakemake script directive
#
# TODO:
# 1) Calibrate against the web interface on a handful of known cases: the
#    predictor was fitted on minfi preprocessNoob M-values, whereas this
#    pipeline produces sesame QCDPB values, so an absolute offset is expected.
#
# History:
#  18-08-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(dplyr))
suppressMessages(library(anndata))
suppressMessages(library(glmnet))
#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input         <- snakemake@input[["adata"]]
    output        <- snakemake@output[["cgc"]]
    predictor_rds <- snakemake@params[["predictor"]]
    predictor_url <- snakemake@params[["predictor_url"]]
    layer         <- snakemake@params[["layer"]]
    value_type    <- snakemake@params[["value_type"]]
}else{
    input         <- '/trinity/home/r115502/SSLOWGRADE/output/Methylation/methylation_data.h5ad'
    output        <- '/trinity/home/r115502/SSLOWGRADE/output/Methylation/CGC/CGC_psi.tsv'
    predictor_rds <- '/trinity/home/r115502/SSLOWGRADE/output/Methylation/CGC/assets/CGC-Psi_predictor_probe_based_lm_v1.0_450k.Rds'
    predictor_url <- 'https://github.com/ErasmusMC-Neuro-Oncology/Continuous_Grading_Classifier/raw/refs/heads/main/assets/CGC-Psi_predictor_probe_based_lm_v1.0_450k.Rds'
    layer         <- NULL
    value_type    <- 'beta'
}
dir.create(dirname(output), recursive = TRUE, showWarnings = FALSE)
#-------------------------------------------------------------------------------
# 0.3 Fetch the predictor
#-------------------------------------------------------------------------------
# The upstream repo ships one glmnet object per platform under assets/. X in
# this pipeline is lifted to HM450 space, so the 450k predictor is the correct
# one even for samples natively run on EPIC/EPICv2.
if(!file.exists(predictor_rds)){
    dir.create(dirname(predictor_rds), recursive = TRUE, showWarnings = FALSE)
    cat(sprintf("[CGC] downloading predictor to %s\n", predictor_rds))
    ok <- tryCatch({
        download.file(predictor_url, predictor_rds, mode = "wb", quiet = TRUE)
        TRUE
    }, error = function(e) FALSE)
    if(!ok || !file.exists(predictor_rds)){
        stop(sprintf(paste0("Could not download the CGC predictor.\n",
                            "On an offline node, fetch it once on the login node:\n",
                            "  curl -L -o %s \\\n    %s"),
                     predictor_rds, predictor_url))
    }
}
predictor <- readRDS(predictor_rds)
#-------------------------------------------------------------------------------
# 1.1 Read data
#-------------------------------------------------------------------------------
adata <- read_h5ad(input)

mat <- if(is.null(layer)) adata$X else adata$layers[[layer]]
mat <- as.matrix(mat)                      # samples x probes
colnames(mat) <- rownames(adata$var)
rownames(mat) <- rownames(adata$obs)

cat(sprintf("[CGC] %d samples x %d probes read from %s\n",
            nrow(mat), ncol(mat), basename(input)))
#-------------------------------------------------------------------------------
# 1.2 Convert to M-values
#-------------------------------------------------------------------------------
# The predictor consumes M-values (minfi ratioConvert what = "M"). X holds beta,
# so convert here; beta is clamped away from 0/1 to keep the logit finite.
if(identical(value_type, 'beta')){
    rng <- range(mat, na.rm = TRUE)
    if(rng[1] < -0.01 || rng[2] > 1.01){
        stop(sprintf(paste0("value_type='beta' but observed range is %.2f..%.2f. ",
                            "Set value_type='M' if this layer already holds M-values."),
                     rng[1], rng[2]))
    }
    eps <- 1e-6
    mat <- log2(pmin(pmax(mat, eps), 1 - eps) / (1 - pmin(pmax(mat, eps), 1 - eps)))
    cat("[CGC] converted beta to M-values\n")
}
#-------------------------------------------------------------------------------
# 1.3 Align to the predictor feature space
#-------------------------------------------------------------------------------
wanted <- rownames(predictor$beta)

# Probes with a zero coefficient contribute nothing, so missingness only
# actually matters for the non-zero set; both are reported separately.
coefs    <- as.matrix(predictor$beta)[, 1]
nonzero  <- names(coefs)[coefs != 0]
missing  <- setdiff(wanted, colnames(mat))

cat(sprintf("[CGC] predictor uses %d probes (%d with non-zero coefficient)\n",
            length(wanted), length(nonzero)))
cat(sprintf("[CGC] %d predictor probes absent from the data, of which %d non-zero\n",
            length(missing), length(intersect(missing, nonzero))))

if(length(intersect(missing, nonzero)) > 0.05 * length(nonzero)){
    warning(sprintf(paste0("%.1f%% of informative probes are absent. CGCpsi ",
                           "values will be biased; check that the h5ad is in ",
                           "HM450 probe space."),
                    100 * length(intersect(missing, nonzero)) / length(nonzero)))
}

data <- matrix(NA_real_, nrow = nrow(mat), ncol = length(wanted),
               dimnames = list(rownames(mat), wanted))
shared <- intersect(wanted, colnames(mat))
data[, shared] <- mat[, shared, drop = FALSE]
#-------------------------------------------------------------------------------
# 1.4 Impute remaining gaps
#-------------------------------------------------------------------------------
# glmnet propagates NA, so pOOBAH-masked probes and liftover gaps have to be
# filled. Per-probe cohort mean where anything was observed; 0 (i.e. beta 0.5)
# for probes with no observation at all, which is arbitrary but bounded - hence
# the per-sample counts written alongside the score.
n_na_before <- rowSums(is.na(data))

probe_means <- colMeans(data, na.rm = TRUE)
probe_means[is.na(probe_means)] <- 0
na_idx <- which(is.na(data), arr.ind = TRUE)
if(nrow(na_idx) > 0){
    data[na_idx] <- probe_means[na_idx[, "col"]]
}

cat(sprintf("[CGC] imputed %d cells (%.2f%% of the predictor matrix)\n",
            sum(n_na_before), 100 * sum(n_na_before) / length(data)))
#-------------------------------------------------------------------------------
# 2.1 Apply the predictor
#-------------------------------------------------------------------------------
pred <- glmnet::predict.glmnet(predictor, data)

if(ncol(pred) > 1){
    warning(sprintf(paste0("Predictor returned %d lambda columns; taking the ",
                           "first. Check whether a specific s= is intended."),
                    ncol(pred)))
}

result <- data.frame(
    sample          = rownames(mat),
    CGC_psi         = as.numeric(pred[, 1]),
    n_probes_used   = length(wanted) - n_na_before,
    n_probes_imputed = as.integer(n_na_before),
    frac_imputed    = round(n_na_before / length(wanted), 4),
    row.names       = NULL,
    check.names     = FALSE
)
#-------------------------------------------------------------------------------
# 3.1 Write output
#-------------------------------------------------------------------------------
write.table(result, output, sep = "\t", quote = FALSE, row.names = FALSE)

cat(sprintf("[CGC] done: CGCpsi range %.3f..%.3f across %d samples\n",
            min(result$CGC_psi), max(result$CGC_psi), nrow(result)))

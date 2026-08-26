#!/usr/bin/env Rscript
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Embeddings_methylation.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# PCA / t-SNE embeddings and all downstream figures, read from the harmonized
# .h5ad written by Preprocess_methylation.R.
#
# Sections 13-14 and the downstream analyses of the original
# Methylation_analysis_new.R. Function internals are verbatim; what changed is
# the ORDER (several objects were used before they were defined) and the removal
# of duplicated definitions.
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv: envs/methylation_embeddings.yaml
# Usage: snakemake script directive
#
# TODO:
# 1)
#
# History:
#  21-08-2026: Split out of Methylation_analysis_new.R
#  22-08-2026: Added purification / LUMP / PC-correlation figures; fixed
#              use-before-definition ordering; removed duplicate definitions
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
library(Rtsne)
library(uwot)
library(RColorBrewer)
library(ggplot2)
library(stringr)
library(dplyr)
library(matrixStats)
library(data.table)
library(anndata)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(cowplot)
library(tibble)
ncores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", 16))
## reticulate >= 1.41 otherwise provisions its own ephemeral uv-managed Python
## env, ignoring the conda env and its anndata/pandas pins.
reticulate::use_condaenv(Sys.getenv("CONDA_PREFIX"), required = TRUE)
#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if (exists("snakemake")) {
  adata_in           <- snakemake@input[['adata']]
  OUT_DIR            <- snakemake@params[['out_dir']]
  embeddings_all_out <- snakemake@output[['embeddings_all']]
  embeddings_sel_out <- snakemake@output[['embeddings_selected']]
  tsne_family_out    <- snakemake@output[['tsne_family']]
  tsne_zoom_out      <- snakemake@output[['tsne_zoom']]
} else {
  adata_in           <- '~/mnt/BIGR_home/MINT/output/Methylation/harmonized/methylation_harmonized.h5ad'
  OUT_DIR            <- '~/mnt/BIGR_home/MINT/output/Methylation/harmonized'
  embeddings_all_out <- 'Embeddings_and_subtypes.csv'
  embeddings_sel_out <- 'Embeddings_and_subtypes_selected.csv'
  tsne_family_out    <- 'tSNE_AllGliomas_w_batch_correction.pdf'
  tsne_zoom_out      <- 'tSNE_AllGliomas_w_batch_correction_zoom.pdf'
}
## Secondary figures are not declared as rule outputs - they land here.
FIG_DIR <- file.path(OUT_DIR, "figures")
LOG_DIR <- file.path(OUT_DIR, "logs")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)
dir.create(LOG_DIR, recursive = TRUE, showWarnings = FALSE)
fig <- function(...) file.path(FIG_DIR, ...)

REF_URL <- "https://raw.githubusercontent.com/livjuergensen/InSilicoPurification/main/data/ref.csv"
RATIOS_TO_USE <- c(0, 0.05, 0.1, 0.15, 0.2, 0.25)
CELLTYPE_SELECTION <- c('T cells', 'Microglia', 'Monocytes')
#-------------------------------------------------------------------------------
# 0.3 Read the harmonized data
#-------------------------------------------------------------------------------
## Rebuilds the probes x samples matrices the analysis expects, plus
## samplesheet / kept_samples / cohort_vector. Matrices come back transposed
## from AnnData's samples x probes layout.
adata <- read_h5ad(adata_in)

samplesheet   <- as.data.frame(adata$obs)
kept_samples  <- rownames(samplesheet)
samplesheet$Sample_ID <- kept_samples
cohort_vector <- samplesheet$Cohort

## ComBat dropped zero-variance probes, stored as NA in the layer and flagged in
## var$in_combat; restrict to those rows so the matrix matches the original.
combat_probes <- rownames(adata$var)[adata$var$in_combat]

mvalues_combat <- t(adata$layers[['mvalue_combat']])
rownames(mvalues_combat) <- rownames(adata$var)
colnames(mvalues_combat) <- rownames(adata$obs)
mvalues_combat <- mvalues_combat[combat_probes, , drop = FALSE]
stopifnot(!anyNA(mvalues_combat))
#-------------------------------------------------------------------------------
# 0.4 Helper functions
#-------------------------------------------------------------------------------
beta_to_mvalue <- function(beta, offset = 1e-6) {
  beta_clipped <- pmin(pmax(beta, offset), 1 - offset)
  log2(beta_clipped / (1 - beta_clipped))
}
mvalue_to_beta <- function(mvalue) {
  b <- 2^mvalue / (2^mvalue + 1)
  pmin(pmax(b, 0), 1)
}
rowVars <- matrixStats::rowVars

## Computed ONCE here. The original derived it late, after several sections had
## already referenced it - LUMP scoring and the Sturm reference signatures both
## need it before the purification block where it used to be defined.
beta_values_combat <- mvalue_to_beta(mvalues_combat)

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

Run_PCA <- function(input_data, n_hvg = 50000){
    beta_matrix_combat <- mvalue_to_beta(input_data)
    probe_vars <- rowVars(beta_matrix_combat)
    if(is.null(n_hvg)){
        N_TOP_VARIABLE_PROBES <- find_elbow(probe_vars)}
    else{
        N_TOP_VARIABLE_PROBES <- as.integer(n_hvg)
        }
    print(paste0('Number of variable probes: ',N_TOP_VARIABLE_PROBES))
    names(probe_vars) <- rownames(beta_matrix_combat)
    top_var_probes <- names(sort(probe_vars, decreasing = TRUE))[
        seq_len(min(N_TOP_VARIABLE_PROBES, length(probe_vars)))
    ]
    mat_top <- beta_matrix_combat[top_var_probes, , drop = FALSE]
    pca_full <- prcomp(t(mat_top), center = TRUE, scale. = FALSE)
    return(pca_full)
}

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

    ## All sample annotation (Methylation_Class, family, Subtype_Pathology, Loc,
    ## Horvath age, NF1_status) is precomputed in obs by
    ## Preprocess_methylation.R, so this is a lookup rather than a rebuild.
    ## Mean_methylation stays here because it is a property of whatever matrix
    ## subset was passed in.
    annot <- samplesheet[match(rownames(pca_scores), samplesheet$Sample_ID), , drop = FALSE]

    embedding_df <- data.frame(
        Sample_ID = rownames(pca_scores),
        PC1 = pca_scores[,1],
        PC2 = pca_scores[,2],
        PC3 = pca_scores[,3],
        PC4 = pca_scores[,4],
        Mean_methylation = colMeans(input_data),
        tSNE1 = tsne_coords[, 1], tSNE2 = tsne_coords[, 2],
        row.names = NULL
    ) %>%
        cbind(annot %>% select(-Sample_ID))

    return(embedding_df)
}

## ---- In-silico purification (subtraction) --------------------------------
## purified = (observed - f*ref)/(1-f), projected back onto the SAME PCA space.
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
      "[%s/%s] Only %d/%d PCA probes covered by the reference - purification may be unreliable for this cell type.",
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

## ---- LUMP purity (Aran, Sirota & Butte 2015, Nat Commun 6:8971) ----------
## purity = mean(beta at 44 leukocyte-unmethylated CpGs) / 0.85
## The pasted IDs had leading zeros stripped (the classic Excel/PDF corruption),
## so every numeric suffix is zero-padded to the standard 8 digits below.
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
lump_cpgs <- sprintf("cg%08d", as.integer(sub("^cg", "", lump_cpgs_raw)))

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

  # Dot radius scales LINEARLY FROM ZERO up to max_r - no minimum floor, so
  # a dot's size is true to its significance: p = 1 (neglogp = 0) renders
  # as a genuine zero-size (invisible) dot, not floored to a visible
  # minimum. Tradeoff versus the old min_r floor: non-significant
  # correlations now disappear from the plot entirely instead of showing
  # as a small dot.
  max_r <- 0.48
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
    matrix(max_r, nrow(neglogp), ncol(neglogp), dimnames = dimnames(neglogp))
  } else {
    (neglogp / max_neglogp_cap) * max_r
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
    rep(max_r, length(legend_breaks))
  } else {
    (legend_breaks / max_neglogp_cap) * max_r
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
    ct <- cor.test(pc_vector[complete_idx], indicator[complete_idx], method ="pearson")
    data.frame(row_label = sprintf("%s = %s", var_name, lvl),
               estimate = unname(ct$estimate), p.value = ct$p.value)
  }))
}
#-------------------------------------------------------------------------------
# 0.5 Colour palettes
#-------------------------------------------------------------------------------
color_pal <- c("#4E79A7","#A0CBE8","#F28E2B","#FFBE7D","#59A14F","#8CD17D","#B6992D","#F1CE63","#499894","#86BCB6","#E15759","#FF9D9A","#79706E","#BAB0AC","#D37295","#FABFD2","#B07AA1","#D4A6C8","#9D7660","#D7B5A6")
names(color_pal) <- c("Blue","Light Blue" ,"Orange","Light Orange","Green","Light Green","Yellow-Green","Yellow","Teal","Light Teal", "Red","Pink","Dark Gray","Light Gray","Pink","Light Pink","Purple","Light Purple","Brown","Light Orange")
cohort_pal <- list('Capper' = color_pal[['Light Gray']],'Lucas' = color_pal[['Red']],'MINT' = color_pal[['Green']], 'Sturm' = color_pal[['Orange']], 'Sturm NF1' = color_pal[['Orange']], 'Sturm NF1wt' = color_pal[['Light Gray']])
loc_pal <- c("#4E79A7","#F28E2B","#76B7B2","#59A14F","#EDC948","#B07AA1","#FF9DA7","#9C755F","#BAB0AC")
subtype_pal <- list('GBM' = color_pal[['Yellow-Green']], DMG = color_pal[['Brown']], 'PA' = color_pal[['Blue']],
                    'RGNT' = color_pal[['Teal']], 'LGG' = color_pal[['Light Teal']], 'IDHmt' = color_pal[['Light Orange']], 'PXA' = color_pal[['Purple']], 'HGAP' = color_pal[['Yellow']],
                    'PGG' = color_pal[['Pink']] , 'IHG' = color_pal[['Red']],'Control' = color_pal[['Dark Gray']], 'Unclass.' = 'white', 'PA, NF1' = color_pal[['Green']])
## Defined once. The original assigned celltype_pal twice, the first from an
## 8-colour seq_along() over all reference types, the second a 4-colour manual
## map named over CELLTYPE_SELECTION - the second silently won, so only it is
## kept here.
celltype_pal <- setNames(color_pal[c('Red','Green','Orange','Teal')][seq_along(CELLTYPE_SELECTION)],
                         CELLTYPE_SELECTION)

## =============================================================================
## Analysis 1.0: All NF1 tumors together with Capper gliomas - tSNE
## =============================================================================
embeddings_all <- Create_embedding_df(mvalues_combat)

write.csv(embeddings_all %>%
            select(Sample_ID,Cohort,Platform,tSNE1,tSNE2,Loc,family,Subtype_Pathology),
            file = embeddings_all_out,row.names = F, quote=F)

## NF1_status comes from obs, so the Sturm samplesheet no longer has to be
## re-read here. Kept as two concatenated pulls rather than one filter(): column
## order into Rtsne determines the layout, and this preserves the original's
## Sturm-first ordering.
Sample_selection <- c(
    samplesheet %>% filter(Cohort == 'Sturm', NF1_status == 'NF1') %>% pull(Sample_ID),
    samplesheet %>% filter(Cohort %in% c('Capper','Lucas','MINT')) %>% pull(Sample_ID)
)

embeddings_all_selected <- Create_embedding_df(mvalues_combat[,Sample_selection]) %>%
    mutate(Cohort = ifelse(Cohort == 'Sturm','Sturm NF1', as.character(Cohort))) %>%
    arrange(Cohort)

write.csv(embeddings_all_selected %>%
            select(Sample_ID,Cohort,Platform,tSNE1,tSNE2,Loc,family,Subtype_Pathology),
            file = embeddings_sel_out,row.names = F, quote=F)

tsne_family <- embeddings_all_selected %>%
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
pdf(fig(tsne_family_out), height = 6 , width = 8)
print(tsne_family)
dev.off()

tsne_family_zoom <- embeddings_all_selected %>%
    filter(tSNE1 > -10, tSNE1 < 10,tSNE2 > -20, tSNE2 < 0 ) %>%
  ggplot(aes(tSNE1, tSNE2, fill = family, shape = Cohort, color = Cohort)) +
  geom_point(alpha = 0.9, size = 4,stroke = 1) +
    scale_fill_manual(values = subtype_pal) +
    scale_color_manual(values = cohort_pal) +
  scale_shape_manual(values = c(22,24,21,23)) +
  theme_classic(base_size = 13) +
  labs(
       fill = 'DKFZ class') +
  guides(
    fill = guide_legend(override.aes = list(shape = 21)),
    shape = guide_legend(override.aes = list(fill = "white"))
  ) +
    theme(legend.position = 'none')

pdf(tsne_zoom_out, height = 5 , width = 5)
print(tsne_family_zoom)
dev.off()

pdf(fig('Figure2_tSNE_unedited.pdf'), height = 9*0.625 , width = 8*2*0.8)
cowplot::plot_grid(tsne_family, tsne_family_zoom, labels= c('A','B') , nrow = 1, rel_widths = c(1,0.75))
dev.off()

## =============================================================================
## Analysis 2.0: All NF1 PAs - PCA, PC associations, in-silico purification
## =============================================================================
PA_samples <- embeddings_all_selected %>%
    filter(Cohort %in% c('MINT','Lucas','Sturm NF1'),
           family %in% c('PA','PA, NF1','Control') |
           (family == 'Unclass.' & Subtype_Pathology == 'PA') |
           (family == 'Control' & Subtype_Pathology == 'PA')) %>%
    pull(Sample_ID) %>% unique()

pca_full_PA   <- Run_PCA(mvalues_combat[,PA_samples])
var_explained <- (pca_full_PA$sdev^2) / sum(pca_full_PA$sdev^2) * 100

## ---- LUMP purity, joined onto the PA embedding ---------------------------
## The original listed LUMP_Immune_pct in sample_vars but never joined it, so
## the missing_cols check below would have stopped the script.
lump_results <- lump_score(beta_values_combat[, PA_samples, drop = FALSE])
write.csv(lump_results, fig('PA_LUMP_purity_immune_pct.csv'), row.names = FALSE)

embeddings_PA <- Create_embedding_df(mvalues_combat[,PA_samples] , perplexity = 10) %>%
    mutate(Cohort = ifelse(Cohort == 'Sturm','Sturm NF1', as.character(Cohort))) %>%
    left_join(lump_results %>% select(Sample_ID, LUMP_Tumor_pct, LUMP_Immune_pct),
              by = 'Sample_ID') %>%
    arrange(Cohort)

Class_probs <- read.delim('~/mnt/BIGR_home/MINT/data/PA_w_class_probs.tsv')
Class_probs <- cbind(
    data.frame(sample = Class_probs$sample,
               PA_score = Class_probs %>% select(contains('PA_')) %>% rowSums()),
    Class_probs %>% select(contains('CTRL')))

embeddings_PA <- embeddings_PA %>% mutate(pct_failed_detection = frac_failed_detection*100 ) %>%
    left_join(Class_probs)


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

pdf(fig('PCA_PA_gliomas_w_batch_correction_family.pdf'), height = 3, width = 6)
print(PCA_PA)
dev.off()

write.csv(embeddings_PA %>%
          mutate(
              idat_green = gsub('/trinity/home/r115502','/home/jurriaan/mnt/BIGR_home/', idat_green),
              idat_red = gsub('/trinity/home/r115502','/home/jurriaan/mnt/BIGR_home/', idat_red),
              idat_green = gsub('.gz','',idat_green),
              idat_red = gsub('.gz','',idat_red),
          ) %>% 
            select(Sample_ID,Cohort,Platform,PC1,PC2,Loc,family,Subtype_Pathology, idat_green,idat_red),
            file = fig('PA_PCA_and_subtypes.csv'),row.names = F, quote=F)

## =============================================================================
## Analysis 2.1: Reference signatures (published + Sturm-derived)
## =============================================================================
## Both references are built BEFORE anything consumes them. In the original the
## Sturm signatures were defined below the block that already used them in
## ref_combined, and beta_values_combat was defined below its first use.
ref_signatures <- read.csv(REF_URL, row.names = 1, check.names = FALSE)
cat("Reference cell types available:", paste(colnames(ref_signatures), collapse = ", "), "\n")
cat(sprintf("Reference covers %d probes\n", nrow(ref_signatures)))

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

## =============================================================================
## Analysis 2.2: In-silico purification of the unclassified PA samples
## =============================================================================
tumor_samples <- embeddings_PA %>% filter(family == 'Unclass.') %>% pull(Sample_ID) %>% unique()
cat(sprintf("Purifying %d unclassified samples\n", length(tumor_samples)))

purification_results_published <- do.call(rbind, lapply(colnames(ref_signatures), function(ct) {
  do.call(rbind, lapply(tumor_samples, function(s) {
    make_purification_series(
      tumor_sample_id = s, ref_signatures = ref_signatures, cell_type = ct,
      beta_matrix = beta_values_combat, pca_object = pca_full_PA,
      fractions = RATIOS_TO_USE
    )
  }))
})) %>% mutate(Reference_Source = "Published (Jurgensen et al.)")

purification_results_sturm <- do.call(rbind, lapply(colnames(ref_signatures_sturm), function(ct) {
  do.call(rbind, lapply(tumor_samples, function(s) {
    make_purification_series(
      tumor_sample_id = s, ref_signatures = ref_signatures_sturm, cell_type = ct,
      beta_matrix = beta_values_combat, pca_object = pca_full_PA,
      fractions = RATIOS_TO_USE
    )
  }))
})) %>% mutate(Reference_Source = "Sturm controls")

purification_results_matched <- rbind(purification_results_published, purification_results_sturm)

## ---- Facet labels, ordered, excluding the 0% start point -------------------
## "-X% normal" rather than "+X% normal" deliberately: this is the opposite
## direction from the additive dilution plot - removing reference signal.
dilution_levels <- paste0('-', sprintf("%d%% normal", round(sort(unique(
  purification_results_matched$fraction_removed[purification_results_matched$fraction_removed > 0]
)) * 100)))
purification_results_matched <- purification_results_matched %>%
  mutate(dilution_label = factor(
    paste0('-', sprintf("%d%% normal", round(fraction_removed * 100))),
    levels = dilution_levels
  ))

arrow_starts_pur <- purification_results_matched %>%
  filter(fraction_removed == 0) %>%
  select(Source_Sample, PC1_start = PC1, PC2_start = PC2) %>%
  unique()

arrow_df_pur <- purification_results_matched %>%
  filter(fraction_removed > 0) %>%
  left_join(arrow_starts_pur, by = "Source_Sample") %>%
  unique()

background_df <- embeddings_PA

p_purification_matched <- ggplot() +
     geom_point(data = background_df %>% filter(family != 'PA'), aes(PC1, PC2),
             color = color_pal[['Dark Gray']], size = 2) +
    geom_point(data = background_df %>% filter(family == 'PA'), aes(PC1, PC2),
               color = color_pal[['Blue']], size = 2) +
  geom_segment(
    data = arrow_df_pur %>% filter(Cell_Type %in% CELLTYPE_SELECTION, PC1_start < 10),
    aes(x = PC1_start, y = PC2_start, xend = PC1, yend = PC2, color = Cell_Type,
        group = interaction(Source_Sample, Cell_Type)),
    arrow = arrow(length = unit(0.12, "cm"), type = "closed"),
    alpha = 0.5, linewidth = 0.25
  ) +
  geom_point(data = arrow_df_pur %>% filter(Cell_Type %in% CELLTYPE_SELECTION, PC1_start < 10, fraction_removed != 0.05),
             aes(PC1, PC2, color = Cell_Type), size = 1, alpha = 0.5) +
  scale_color_manual(values = celltype_pal) +
  facet_wrap(~dilution_label, nrow=1) +
  theme_classic(base_size = 13) +
    labs(
         x = paste0('PC1 (',round(var_explained[1]),'%)' ),
        y = paste0('PC2 (',round(var_explained[2]),'%)' ), color = "Jürgensen et al reference") +
    theme(legend.position = 'bottom')

pdf(fig('Purification_PA_published_and_sturm_references.pdf'), height = 3, width = 8)
print(p_purification_matched)
dev.off()


## =============================================================================
## Analysis 3.1: SAMPLE-LEVEL - PC1/PC2 vs Horvath, Cohort, LUMP, Loc, etc.
## =============================================================================
## Every variable in sample_vars is handled automatically by
## compute_pc_pearson_rows(), regardless of type - no need to pre-process,
## dummy-code, or separate continuous from categorical.
sample_vars <- c( "Horvath", "Cohort", "Platform", "LUMP_Immune_pct", 'Loc','Mean_methylation','pct_failed_detection','PA_score','CTRL_REACTIVE')
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
rownames(cor_mat_sample)[rownames(cor_mat_sample) == 'pct_failed_detection'] <- '%-age failed probes'
rownames(cor_mat_sample)[rownames(cor_mat_sample) == 'CTRL_REACTIVE'] <- 'Cal. score, control reactive TME'
rownames(cor_mat_sample)[rownames(cor_mat_sample) == 'PA_score'] <- 'Cal. score, PA (cumulative)'

rownames(pval_mat_sample) <- rownames(cor_mat_sample)

pdf(fig("PC_correlation_dotheatmap_sample_level.pdf"), height = 3.5, width = 5)
plot_corr_dotplot(cor_mat_sample, pval_mat_sample)
dev.off()

## ---- Save the underlying numbers too, not just the figures -----------------
write.csv(as.data.frame(cor_mat_sample) %>% tibble::rownames_to_column("Variable"),
          fig("PC_correlation_sample_level_r.csv"), row.names = FALSE)
write.csv(as.data.frame(pval_mat_sample) %>% tibble::rownames_to_column("Variable"),
          fig("PC_correlation_sample_level_pval.csv"), row.names = FALSE)

## =============================================================================
## Analysis 4.0: COMBINED FIGURES
## =============================================================================
## plot_corr_dotplot() draws through ComplexHeatmap, which writes straight to
## the device and returns nothing cowplot can use. draw_now = FALSE returns the
## pieces instead, and grid.grabExpr() captures the drawing as a grob.
built <- plot_corr_dotplot(cor_mat_sample, pval_mat_sample, draw_now = FALSE)
heatmap_grob <- grid.grabExpr(
  draw(built$ht, annotation_legend_list = built$legend_list, merge_legend = TRUE)
)

## Two panels: PCA + correlation dotplot
pdf(fig('PCA_and_dotplot.pdf'), width = 10, height = 3)
cowplot::plot_grid(PCA_PA, heatmap_grob, labels = c('A','B'))
dev.off()

## Three panels: PCA + correlation dotplot + purification
pdf(fig('Figure3_PA_PCA_and_purification.pdf'), width = 10, height = 3*2)
cowplot::plot_grid(
  cowplot::plot_grid(PCA_PA, heatmap_grob, labels = c('A','B')),
  p_purification_matched, labels = c('','C'), nrow = 2, rel_heights = c(1,0.9))
dev.off()

cat(sprintf("[%s] Done. Figures written to %s\n", Sys.time(), FIG_DIR))



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


purified_probs <- read.delim('~/mnt/BIGR_home/MINT/data/idat/purified/heidelberg/Epignostix_probabilities.csv', sep = ',')


data.frame(fractions = as.integer(gsub('M06_pur','',purified_probs$sample)),
           PA_score = purified_probs %>% select(contains('PA')) %>% rowSums(),
           Control_score = purified_probs %>% select(contains('CTRL')) %>% rowSums()) %>%
    tidyr::pivot_longer(cols = c(PA_score,Control_score)) %>%
    ggplot(aes(fractions,value,color = name)) + geom_line()




## =============================================================================
## Analysis 5.0: ANALYSIS WITH NON-NF1 PAs
## =============================================================================
PA_samples_non_NF1_Sturm  <- embeddings_all %>% filter(Cohort == 'Sturm', NF1_status == 'NF1wt',Subtype_Pathology == 'PA',family == 'PA') %>% pull(sample)
PA_samples_extended <- c(PA_samples,PA_samples_non_NF1_Sturm)
lump_results <- lump_score(beta_values_combat[, PA_samples_extended, drop = FALSE])

embeddings_PA_extended <- Create_embedding_df(mvalues_combat[,PA_samples_extended]) %>%
    mutate(Cohort = case_when(
               Cohort== 'Sturm' & NF1_status == 'NF1wt' ~ 'Sturm NF1wt',
               Cohort== 'Sturm' & NF1_status == 'NF1' ~ 'Sturm NF1',
               TRUE ~ Cohort)) %>%
    left_join(lump_results %>% select(Sample_ID, LUMP_Tumor_pct, LUMP_Immune_pct),
              by = 'Sample_ID')


pca_full_PA   <- Run_PCA(mvalues_combat[,PA_samples_extended])
var_explained <- (pca_full_PA$sdev^2) / sum(pca_full_PA$sdev^2) * 100

PCA_PA_extended <- embeddings_PA_extended %>%
    arrange(desc(Cohort)) %>%
    ggplot(aes(PC1, PC2, fill = family, shape = Cohort, color = Cohort)) +
    geom_point( size = 2.5,stroke = 0.9) +
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
    geom_point(data = embeddings_PA_extended %>% filter(Methylation_Class == 'PA, NF1-associated'),aes(PC1,PC2), shape = 8, color = 'black', inherit.aes = F, size = 0.9)


pdf(fig('PCA_PA_extended_gliomas_w_batch_correction_family.pdf'), height = 3, width = 6)
print(PCA_PA_extended)
dev.off()

tSNE_PA_extended <- embeddings_PA_extended %>%
    arrange(desc(Cohort)) %>%
    ggplot(aes(tSNE1, tSNE2, fill = family, shape = Cohort, color = Cohort)) +
    geom_point( size = 2.5,stroke = 0.9) +
    scale_fill_manual(values = subtype_pal) +
    scale_color_manual(values = cohort_pal) +
    scale_shape_manual(values = c(22,24,21,23)) +
    theme_classic(base_size = 13) +
    labs(
        fill = 'DKFZ class',
        x = 'tSNE1',
        y= 'tSNE2') + 
    guides(
        fill = guide_legend(override.aes = list(shape = 21))) +
    geom_point(data = embeddings_PA_extended %>% filter(Methylation_Class == 'PA, NF1-associated'),aes(tSNE1,tSNE2), shape = 8, color = 'black', inherit.aes = F, size = 0.9)


pdf(fig('tSNE_PA_extended_gliomas_w_batch_correction_family.pdf'), height = 3, width = 6)
print(tSNE_PA_extended)
dev.off()




library(limma)
mvalues_PA_extended <- mvalues_combat[, PA_samples_extended, drop = FALSE]
lump_covariate <- lump_results$LUMP_Immune_pct[match(colnames(mvalues_PA_extended), lump_results$Sample_ID)]

mvalues_PA_extended_lumpRegressed <- removeBatchEffect(
    mvalues_PA_extended,
    covariates = lump_covariate
)


pca_full_PA   <- Run_PCA(mvalues_PA_extended_lumpRegressed[,PA_samples_extended])
var_explained <- (pca_full_PA$sdev^2) / sum(pca_full_PA$sdev^2) * 100

embeddings_PA_extended_lumpRegressed <- Create_embedding_df(mvalues_PA_extended_lumpRegressed) %>%
    mutate(Cohort = case_when(
               Cohort== 'Sturm' & NF1_status == 'NF1wt' ~ 'Sturm NF1wt',
               Cohort== 'Sturm' & NF1_status == 'NF1' ~ 'Sturm NF1',
               TRUE ~ Cohort)) %>%
    left_join(lump_results %>% select(Sample_ID, LUMP_Tumor_pct, LUMP_Immune_pct),
              by = 'Sample_ID')

add_inside_title <- function(p, label, hjust = -0.05, vjust = 1.6, size = 4, fontface = "bold") {
    p + annotate("text", x = -Inf, y = Inf, label = label,
                 hjust = hjust, vjust = vjust, size = size, fontface = fontface)
}

## ---- Rebuild the two "after" plots WITHOUT labs(title=...) - the panel
## itself is now identical in size/position to the "before" plots, since
## there's no title row competing for space.
PCA_PA_extended_lumpRegressed <- embeddings_PA_extended_lumpRegressed %>%
    arrange(desc(Cohort)) %>%
    ggplot(aes(PC1, PC2, fill = family, shape = Cohort, color = Cohort)) +
    geom_point(size = 2.5, stroke = 0.9) +
    scale_fill_manual(values = subtype_pal) +
    scale_color_manual(values = cohort_pal) +
    scale_shape_manual(values = c(22, 24, 21, 23)) +
    theme_classic(base_size = 13) +
    labs(fill = 'DKFZ class',
         x = paste0('PC1 (', round(var_explained[1]), '%)'),
         y = paste0('PC2 (', round(var_explained[2]), '%)'),
         subtitle = 'Corrected for LUMP %-age') +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    geom_point(data = embeddings_PA_extended_lumpRegressed %>% filter(Methylation_Class == 'PA, NF1-associated'),
               aes(PC1, PC2), shape = 8, color = 'black', inherit.aes = FALSE, size = 0.9)

tSNE_PA_extended_lumpRegressed <- embeddings_PA_extended_lumpRegressed %>%
    arrange(desc(Cohort)) %>%
    ggplot(aes(tSNE1, tSNE2, fill = family, shape = Cohort, color = Cohort)) +
    geom_point(size = 2.5, stroke = 0.9) +
    scale_fill_manual(values = subtype_pal) +
    scale_color_manual(values = cohort_pal) +
    scale_shape_manual(values = c(22, 24, 21, 23)) +
    theme_classic(base_size = 13) +
    labs(fill = 'DKFZ class', x = 'tSNE1', y = 'tSNE2',subtitle = 'Corrected for LUMP %-age') +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    geom_point(data = embeddings_PA_extended_lumpRegressed %>% filter(Methylation_Class == 'PA, NF1-associated'),
               aes(tSNE1, tSNE2), shape = 8, color = 'black', inherit.aes = FALSE, size = 0.9)

## ---- Strip legends AND add the inside-panel titles ------------------------
no_legend <- theme(legend.position = "none")

PCA_PA_extended_nl                <- PCA_PA_extended                + no_legend
tSNE_PA_extended_nl               <- tSNE_PA_extended               + no_legend
PCA_PA_extended_lumpRegressed_nl  <- PCA_PA_extended_lumpRegressed + no_legend
tSNE_PA_extended_lumpRegressed_nl <- tSNE_PA_extended_lumpRegressed + no_legend

## ---- Shared horizontal legend (unchanged from before) ---------------------
shared_legend <- cowplot::get_legend(
    PCA_PA_extended +
    guides(fill = guide_legend(override.aes = list(shape = 21), nrow = 4)))


## ---- 2x2 grid + legend row --------------------------------------------------
grid_2x2 <- cowplot::plot_grid(
    PCA_PA_extended_nl, tSNE_PA_extended_nl,
    PCA_PA_extended_lumpRegressed_nl, tSNE_PA_extended_lumpRegressed_nl,
    nrow = 2, labels = c('A', 'B', 'C', 'D'),
    align = "hv"   # now that no plot has a labs(title), this aligns panels
                    # (not just outer plot edges) across the whole grid
)

final_grid <- cowplot::plot_grid(
     grid_2x2,shared_legend,
    ncol = 2, rel_widths =  c(1,0.2)
)

pdf(fig('Figure4_PCA_tSNE_PA_extended_before_after_LUMP_grid.pdf'), height = 5, width = 8)
print(final_grid)
dev.off()






pca_full_PA_extended_lumpRegressed <- Run_PCA(mvalues_PA_extended)#_lumpRegressed)
var_explained_lumpRegressed <- (pca_full_PA_extended$sdev^2) /
    sum(pca_full_PA_extended_lumpRegressed$sdev^2) * 100

embeddings_PA_extended_lumpRegressed <- data.frame(
    Sample_ID = rownames(pca_full_PA_extended_lumpRegressed$x),
    PC1 = pca_full_PA_extended_lumpRegressed$x[, 1],
    PC2 = pca_full_PA_extended_lumpRegressed$x[, 2]
) %>%
    left_join(
        embeddings_PA_extended %>%
            select(Sample_ID, Cohort, family, Subtype_Pathology, LUMP_Tumor_pct, LUMP_Immune_pct),
        by = "Sample_ID"
    )




plot_data <- cbind(embeddings_PA_extended_lumpRegressed %>% select(Sample_ID,family,Cohort), pca_full_PA_extended_lumpRegressed$x[,1:50],lump_covariate)
pc_cols <- grep("^PC[0-9]+$", colnames(plot_data), value = TRUE) 
pc_matrix <- as.matrix(plot_data[, pc_cols])
rownames(pc_matrix) <- plot_data$Sample_ID



sample_cor <- cor(t(pc_matrix), method = "pearson")  # samples x samples
 
## ---- 4. Annotation data, in the SAME sample order as sample_cor ----------
annotation_df <- plot_data[match(rownames(sample_cor), plot_data$Sample_ID),
                            c("Cohort", "family",'lump_covariate')]




top_anno <- HeatmapAnnotation(
  Cohort = annotation_df$Cohort,
  `DKFZ class` = annotation_df$family,
  `LUMP leukocyte %-age` = lump_covariate,
  annotation_name_gp = gpar(fontsize = 8),
  col = list(Cohort = unlist(cohort_pal), `DKFZ class` = unlist(subtype_pal),
              `LUMP leukocyte %-age` = colorRamp2(c(min(lump_covariate) ,mean(lump_covariate), max(lump_covariate)), c("#4E79A7", "white", "#F28E2B"))),
  which = "column",
  show_legend = TRUE
)



left_anno <- HeatmapAnnotation(
  Cohort = annotation_df$Cohort,
  `DKFZ class` = annotation_df$family,
    `LUMP leukocyte %-age` = lump_covariate,
  annotation_name_gp = gpar(fontsize = 8),
  col = list(Cohort = unlist(cohort_pal), `DKFZ class` = unlist(subtype_pal),
                           `LUMP leukocyte %-age` = colorRamp2(c(min(lump_covariate) ,mean(lump_covariate), max(lump_covariate)), c("#4E79A7", "white", "#F28E2B"))),
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
ht_sample_corr 




# Side-by-side comparison, same plotting convention as PCA_PA_extended
PCA_PA_extended_lumpRegressed <- embeddings_PA_extended_lumpRegressed %>%
    arrange(desc(Cohort)) %>%
    ggplot(aes(PC1, PC2, fill = family, shape = Cohort, color = Cohort)) +
    geom_point(size = 3, stroke = 1) +
    scale_fill_manual(values = subtype_pal) +
    scale_color_manual(values = cohort_pal) +
    scale_shape_manual(values = rev(c(22, 24, 21, 23))) +
    theme_classic(base_size = 13) +
    labs(
        fill = "DKFZ class",
        x = paste0("PC1 (", round(var_explained_lumpRegressed[1]), "%)"),
        y = paste0("PC2 (", round(var_explained_lumpRegressed[2]), "%)"),
        title = "After regressing out LUMP_Immune_pct"
    ) +
    guides(fill = guide_legend(override.aes = list(shape = 21)))


pdf("PCA_PA_extended_before_after_LUMP_regression.pdf", height = 3.5, width = 12)
cowplot::plot_grid(
    PCA_PA_extended + labs(title = "Before"),
    PCA_PA_extended_lumpRegressed,
    labels = c("A", "B")
)
dev.off()

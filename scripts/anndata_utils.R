#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# anndata_utils.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Utils functions for working with anndata in R
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv: R
# Usage:
#
# adata<- read_h5ad('~/Projects/MINT/MethylationArray-snake/output/methylation/methylation_data.h5ad')
# adata <- Identify_hv_sites(adata)
# adata <- Run_PCA(adata, layer = 'beta')
# adata <- Run_TSNE(adata, perplexity = 10)
#
# TODO:
# 1) 
#
# History:
#  12-03-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(anndata))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(ComplexHeatmap))

#-------------------------------------------------------------------------------
# 1 Modify anndata
#-------------------------------------------------------------------------------
# Identify highly variable cpg sites
Identify_hv_sites <- function(adata, layer = 'X',top_fraction = 0.2){
    mat <- get_matrix(adata, layer)
    # calculate variance
    var <- apply(mat, 2, var, na.rm = TRUE)
    names(var) <- adata$var_names
    # Select top fraction
    n_top <- floor(top_fraction * ncol(mat)) 
    highly_variable_sites <- names(sort(var, decreasing = TRUE))[1:n_top]
    # add column
    adata$var$is_highly_variable <- adata$var_names %in% highly_variable_sites
    return(adata)
}


# Run PCA 
Run_PCA <- function(adata,layer = 'X',scaling = T,use_highly_variable=T){
    mat <- get_matrix(adata, layer)
    # subset to highly variable 
    if(use_highly_variable == T){
        mat <- mat[,adata$var$is_highly_variable]
    }
    # scale data
    if(scaling == T){
        mat <- scale(mat)
    }
    # compute pca model
    pca_mod <- prcomp(mat, center = F)

    # extract scores/loadings and save in adata.uns
    adata$uns[['PCA_scores']] = as.data.frame(pca_mod$x) 
    adata$uns[['PCA_loadings']] = as.data.frame(pca_mod$rotation )

    return(adata)
}

# Run TSNE 
Run_TSNE <- function(adata,perplexity, layer = 'PCA_scores', use_correlations=T){
    mat <- get_matrix(adata, layer)

    if(use_correlations == T){
        mat <- cor(t(mat))
        }
    # calculate 
    # compute tsne model
    tsne_mod <- Rtsne::Rtsne(mat, perplexity = perplexity)

    # extract coordinates
    tsne_coordinates <- as.data.frame(tsne_mod$Y)
    rownames(tsne_coordinates) <- adata$obs_names
    colnames(tsne_coordinates) <- c('tSNE1','tSNE2')
    
    # extract scores/loadings and save in adata.uns
    adata$uns[['tSNE_coords']] = tsne_coordinates

    return(adata)
}


plot(adata$X[1,],
     adata$uns[['beta']][1,]
     )


adata

adata$layers[['beta']] %>% dim()


plot(adata$X[2,], adata$layers[['beta']][2,])

#-------------------------------------------------------------------------------
# 2  Plot adata
#-------------------------------------------------------------------------------
function(adata, layer = 'X'){}


#-------------------------------------------------------------------------------
# 3  Misc functions
#-------------------------------------------------------------------------------
# Helper function to extract data matrix from adata
get_matrix <- function(adata, layer) {
    if (layer == "X") {
        mat <- adata$X
        rownames(mat) <- adata$obs$sample
    } else if (layer %in% names(adata$layers)) {
        mat <- adata$layers[[layer]]
    } else if (layer %in% names(adata$uns)) {
        mat <- adata$uns[[layer]]
    } else {
        stop(sprintf("Layer '%s' not found in adata$X, adata$layers, or adata$uns", layer))
    }
    return(mat)
}

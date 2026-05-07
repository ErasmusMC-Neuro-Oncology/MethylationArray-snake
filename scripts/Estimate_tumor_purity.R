#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Estimate_tumor_purity.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Run InfiniumPurify and RF_purtiy on methylation data
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv: methylation
# Usage:
#
#
# TODO:
# 1) 
#
# History:
#  18-03-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Setup and load packages
#-------------------------------------------------------------------------------
source('scripts/install_github_packages.R')
suppressMessages(library(anndata))
suppressMessages(library(dplyr))
suppressMessages(library(RFpurify))
suppressMessages(library(InfiniumPurify))

source('scripts/anndata_utils.R')
# Use correct reticulate environment
reticulate::use_condaenv(Sys.getenv("CONDA_PREFIX"), required = TRUE)
#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input<- snakemake@input[[1]]
    tumor_type <- snakemake@params[['tumor_type']]
    output <- snakemake@output[[1]]
}else{
    input <- '/home/jurriaan/Projects/Capper_Methylation/MethylationArray-snake/output/methylation/methylation_data_MINT.h5ad'
    output <- '/home/jurriaan/Projects/Capper_Methylation/MethylationArray-snake/output/results/Tumor_purities.txt'
}
#-------------------------------------------------------------------------------
# 1.1 Read data
#-------------------------------------------------------------------------------
# Read adata
adata <- read_h5ad(input)

#-------------------------------------------------------------------------------
# 2.1 Run InfiniumPurify
#-------------------------------------------------------------------------------
# fetch beta matrix
beta <- get_matrix(adata, 'beta')
# Run InfiniumPurify
InfiniumPurify_purity <- data.frame(
    sample = adata$obs$sample,
    InfiniumPurify_purity_LGG= InfiniumPurify::getPurity(t(beta), tumor.type = 'LGG'),
    InfiniumPurify_purity_GBM= InfiniumPurify::getPurity(t(beta), tumor.type = 'GBM'))

#-------------------------------------------------------------------------------
# 2.2 Run RFpurify
#-------------------------------------------------------------------------------
# fetch beta matrix
beta <- get_matrix(adata, 'beta_raw')
# fetch feature matrix and impute missing values
featuremat_ABSOLUTE <- beta[match(rownames(RFpurify_ABSOLUTE$importance), rownames(beta)), , drop = FALSE]
featuremat_ESTIMATE <- beta[match(rownames(RFpurify_ESTIMATE$importance), rownames(beta)), , drop = FALSE]
featuremat_ABSOLUTE[is.na(featuremat_ABSOLUTE)] <- 0.5
featuremat_ESTIMATE[is.na(featuremat_ESTIMATE)] <- 0.5

RFpurify_purity <- data.frame(
    sample = adata$obs$sample,
    RFpurity_ABSOLUTE = predict(RFpurify_ABSOLUTE, t(featuremat_ABSOLUTE)),
    RFpurity_ESTIMATE = predict(RFpurify_ESTIMATE, t(featuremat_ESTIMATE)))

#-------------------------------------------------------------------------------
# 3.1 Join an write to file
#-------------------------------------------------------------------------------
InfiniumPurify_purity %>%
    left_join(RFpurify_purity) %>%
    write.table(output, sep = '\t', quote = F, row.names = F)

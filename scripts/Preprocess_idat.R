#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Preprocess_idat.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Preprocess .idat files
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
#  11-03-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(dplyr))
suppressMessages(library(minfi))

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input<- snakemake@input
    output <- snakemake@output
}else{
    input <-  '../MINT/data/samplesheets/samplesheet_methylation.csv'
    output <- 'output/'
}

#-------------------------------------------------------------------------------
# 1.1 Read data
#-------------------------------------------------------------------------------
# Read samplesheet
samplesheet <- read.delim(input , sep = ',')  %>% mutate(idat_basename = gsub("_Red.idat$", "", idat_red))

# Read idats
raw_intensity_data <- read.metharray(samplesheet$idat_basename, force=T)
# Add sample IDs
colnames(raw_intensity_data) <- samplesheet$sample
#-------------------------------------------------------------------------------
# 2.1 Quality Control: filter out samples with >5% failed probes
#-------------------------------------------------------------------------------


#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create_SampleData.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Combine MethylationArray pipeline results into sample data
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv: 
# Usage: 
#
# TODO:
#
# History:
#  20-05-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(dplyr))

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_Classes <- snakemake@input[["Classes"]]
    input_purities <- snakemake@input[["purities"]]
    output <- snakemake@output[[1]]
    
}else{
    input_Classes <- '/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/Methylation/results/Methylation_Classes.txt'
    input_purities <-'/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/Methylation/results/Tumor_purities.txt'
    output <- "/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/Methylation/sampledata/SampleData_Methylation.txt"
}



#-------------------------------------------------------------------------------
# 1.1 Read data 
#-------------------------------------------------------------------------------
# Read datasets
classes <- read.delim(input_Classes)
purities <- read.delim(input_purities)
#-------------------------------------------------------------------------------
# 1.2 Reformat and join data
#-------------------------------------------------------------------------------
# Join data
SampleData <- classes %>%
    left_join(purities %>% mutate(sample = gsub('_tumor1','',sample)))
   
#-------------------------------------------------------------------------------
# 2.1 Write to file
#-------------------------------------------------------------------------------
write.table(SampleData, file =  output, sep = '\t', quote = F, row.names = F)

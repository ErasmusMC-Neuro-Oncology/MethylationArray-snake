#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Create_SampleData.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Fetch Methylation array results
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
#  18-08-2026: Compile results with new adata
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(dplyr))

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_predictions <- snakemake@input[["predictions"]]
    input_scores <- snakemake@input[["scores"]]
    input_cgc <- snakemake@input[["cgc"]]
    output <- snakemake@output[[1]]    
}else{
    input_predictions <- '/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/Methylation/crossNN/crossNN_predictions.tsv'
    input_scores <- '/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/Methylation/crossNN/crossNN_scores.tsv'
    input_cgc <- '/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/Methylation/CGC/CGC_psi.tsv'
    output <- "/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/Methylation/sampledata/SampleData_Methylation.txt"
}



#-------------------------------------------------------------------------------
# 1.1 Read data 
#-------------------------------------------------------------------------------
# Read datasets
predictions <- read.delim(input_predictions)
scores <- read.delim(input_scores)
cgc <- read.delim(input_cgc)
#-------------------------------------------------------------------------------
# 1.2 Reformat and join data
#-------------------------------------------------------------------------------
# Join data
SampleData <- cgc %>%
    left_join(predictions) %>%
    mutate(CGC = log(scores$A.IDH..HG / scores$A.IDH),
           CGC_curated = case_when(
               grepl('O IDH',predicted_class) ~ CGC_psi,
               grepl('A IDH',predicted_class) ~ CGC,
               TRUE ~ NA)) %>%
    mutate(patient = gsub('_R1_tumor1|_tumor1','',sample))


delta_CGC <- SampleData %>%
    mutate(Surgery = ifelse(grepl('R1',sample),'primary','recurrence')) %>%
    tidyr::pivot_wider(id_cols = patient,names_from = Surgery,values_from = CGC_curated) %>%
    mutate(delta_CGC = recurrence-primary) %>% select(patient,delta_CGC)

SampleData <- SampleData %>% left_join(delta_CGC) %>%
    select(sample,predicted_class,score,interpretation,CGC_psi,CGC,CGC_curated,delta_CGC)

#-------------------------------------------------------------------------------
# 2.1 Write to file
#-------------------------------------------------------------------------------
write.table(SampleData, file =  output, sep = '\t', quote = F, row.names = F)

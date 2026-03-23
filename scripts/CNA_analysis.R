#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# CNA_analysis.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Perform CNA analysis from methylation arrays
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv:  methylation
# Usage: 
#
# TODO:
# 1) 
#
# History:
#  17-03-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
if(!"IlluminaHumanMethylationEPICmanifest" %in% installed.packages()){devtools::install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38")}
if(!"IlluminaHumanMethylationEPICanno.ilm10b5.hg38" %in% installed.packages()){devtools::install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38")}
if(!'conumee2' %in% installed.packages()){ devtools::install_github("hovestadtlab/conumee2", subdir = "conumee2")}
suppressMessages(library(dplyr))
suppressMessages(library(minfi))
suppressMessages(library(anndata))
suppressMessages(library(conumee2))

# Use correct reticulate environment
reticulate::use_condaenv(Sys.getenv("CONDA_PREFIX"), required = TRUE)

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_query <- snakemake@input[['query']]
    input_reference <- snakemake@input[['reference']]
    output_segmented <-  snakemake@output[['Segmented']]
    profile_dir  <-  snakemake@output[['Profile_dir']]
}else{
    input_query <- '/home/jurriaan/Projects/Capper_Methylation/MethylationArray-snake/output/methylation/MINT/methylation_object.Rds'
    input_reference <- '/home/jurriaan/Projects/Capper_Methylation/MethylationArray-snake/output/methylation/Pai/methylation_object.Rds'
    output_segmented <- '/home/jurriaan/Projects/Capper_Methylation/MethylationArray-snake/output/CNAs/MINT/Segmented_CNAs.txt'
    profile_dir <- '/home/jurriaan/Projects/Capper_Methylation/MethylationArray-snake/output/CNAs/MINT/plots/'
}
#-------------------------------------------------------------------------------
# 1.1 Read data
#-------------------------------------------------------------------------------
# Read methylation data
query <- readRDS(input_query)
reference <- readRDS(input_reference)
# Combine arrays
methylation_data <- minfi::combineArrays(query, reference)

array <- strsplit(annotation(methylation_data)[1],'IlluminaHumanMethylation')[[1]][2]

#-------------------------------------------------------------------------------
# 3.1 CNV analysis
#-------------------------------------------------------------------------------
# Fetch annotations and subset probes
data(exclude_regions)
data(detail_regions)
anno <- CNV.create_anno(array_type = array, exclude_regions = exclude_regions, detail_regions =detail_regions)
anno@probes <- anno@probes[names(anno@probes) %in% rownames(methylation_data)]



#-------------------------------------------------------------------------------
# create CNV object
CNV_object <- CNV.load(methylation_data[,1:ncol(query)])
CNV_control <- CNV.load(methylation_data[,(ncol(query)+1):ncol(methylation_data)])

# Estimate CNVs 
CNVs <- CNV.fit(CNV_object, CNV_control , anno)

# Perform binning and segmentation
CNVs <- CNV.bin(CNVs)
CNVs <- CNV.segment(CNVs)



#-------------------------------------------------------------------------------
# 3.2 Plot profiles
#-------------------------------------------------------------------------------

CNV.genomeplot(CNVs[1])




#-------------------------------------------------------------------------------
# 3.3 Create CNV export
#-------------------------------------------------------------------------------

CNV.write(CNVs, what = "segments", file = output_segmented )


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
suppressMessages(library(dplyr))
suppressMessages(library(minfi))
suppressMessages(library(anndata))
suppressMessages(library(conumee))

# Use correct reticulate environment
reticulate::use_condaenv(Sys.getenv("CONDA_PREFIX"), required = TRUE)

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input <- snakemake@input[[1]]
    output <-  snakemake@output[['Segmented']]
    profile_dir  <-  snakemake@output[['Profile_dir']]
}else{
    input <- '/home/jurriaan/Projects/MINT/data/samplesheets/samplesheet_methylation.csv'
    output <- '/home/jurriaan/Projects/Capper_Methylation/MethylationArray-snake/output/CNAs/Segmented_CNAs_MINT.txt'
    profile_dir <- '/home/jurriaan/Projects/Capper_Methylation/MethylationArray-snake/output/CNAs/plots/MINT/' 
}
#-------------------------------------------------------------------------------
# 1.1 Read data
#-------------------------------------------------------------------------------
# Read samplesheets
samplesheet_query <- read.delim(input , sep = ',')  %>%
    mutate(idat_basename = gsub("_Red.idat$", "", idat_red),
           group = 'query') 

# Fetch Pai et al non-tumor samples
samplesheet_reference <-
    data.frame(idat_red = list.files('/data/Resources/datasets/Pai/idat/',pattern = 'Red', full.names = T)) %>%
    mutate(idat_basename = gsub("_Red.idat$", "", idat_red),
           sample = basename(idat_basename),
           group = 'reference')

# Combine samplesheets
samplesheet <- rbind(
    samplesheet_query %>% select(sample,idat_basename,group),
    samplesheet_reference %>% select(sample,idat_basename,group))

# Read idats
raw_intensity_data <- read.metharray(samplesheet$idat_basename, force=T, verbose = T)

#-------------------------------------------------------------------------------
# 2.2 Normalization: Perform Noob normalization
#-------------------------------------------------------------------------------
colnames(raw_intensity_data) <- samplesheet$sample
normalized_data <- preprocessNoob(raw_intensity_data)

# Keep only probes that succeeded in all samples
detection_pvalues <-  detectionP(raw_intensity_data)
keep_probes <- rowSums(detection_pvalues < 0.01) == ncol(normalized_data)
normalized_data <- normalized_data[keep_probes, ]

array <- strsplit(annotation(raw_intensity_data)[1],'IlluminaHumanMethylation')[[1]][2]
#-------------------------------------------------------------------------------
# 3.1 CNV analysis
#-------------------------------------------------------------------------------
# Fetch annotations and subset probes
data(exclude_regions)
data(detail_regions)
anno <- CNV.create_anno(array_type = array, exclude_regions = exclude_regions, detail_regions =detail_regions)
anno@probes <- anno@probes[names(anno@probes) %in% rownames(normalized_data)]

# create CNV object
CNV_object <- CNV.load(normalized_data)

# Estimate CNVs and save objects in list
tumor_samples <- which(samplesheet$group == 'query')
cnv_list <- lapply(tumor_samples, function(i) {
    CNV.fit(CNV_object[i,],CNV_object[samplesheet$group == 'reference',],anno=anno )
})

# Perform binning and segmentation
cnv_list <- lapply(cnv_list, CNV.bin)
cnv_list <- lapply(cnv_list, CNV.segment)

#-------------------------------------------------------------------------------
# 3.2 Create CNV export
#-------------------------------------------------------------------------------



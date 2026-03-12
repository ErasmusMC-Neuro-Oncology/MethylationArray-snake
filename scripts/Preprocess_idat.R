#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Preprocess_idat.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Preprocess .idat files
#
# Author: Jurriaan Janssen (j.janssen.1@erasmusmc.nl)
#
# condaenv: /home/jurriaan/Projects/MethylationArray-snake/.snakemake/conda/f85d34e2
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
if(!"IlluminaHumanMethylationEPICmanifest" %in% installed.packages()){devtools::install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38")}
if(!"IlluminaHumanMethylationEPICanno.ilm10b5.hg38" %in% installed.packages()){devtools::install_github("achilleasNP/IlluminaHumanMethylationEPICanno.ilm10b5.hg38")}
suppressMessages(library(dplyr))
suppressMessages(library(minfi))
suppressMessages(library(anndata))

# Use correct reticulate environment
reticulate::use_condaenv(Sys.getenv("CONDA_PREFIX"), required = TRUE)

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input<- snakemake@input[[1]]
    Zhou_input <- snakemake@params[['Zhou_probes']]
    CrossReactive_input <- snakemake@params[['CrossReactive_probes']]
    Problematic_input <- snakemake@params[['Problematic_probes']]
    output <- snakemake@output[[1]]
}else{
    input <-  '../MINT/data/samplesheets/samplesheet_methylation.csv'
    Zhou_input <- '/data/Resources/EPIC/manifest/AppendixD_Zhou_et_al_MASKgeneral_list.txt'
    CrossReactive_input <- '/data/Resources/EPIC/manifest/AppendixE_CrossReactiveProbes_EPICv1.txt'
    Problematic_input <- '/data/Resources/EPIC/manifest/AppendixF_ProblematicProbes_EPICv1-b5.txt'
    output <- 'output/methylation/methylation_data.h5ad'
}
#-------------------------------------------------------------------------------
# 1.1 Read data
#-------------------------------------------------------------------------------
# Read samplesheet
samplesheet <- read.delim(input , sep = ',')  %>%
    mutate(idat_basename = gsub("_Red.idat$", "", idat_red),
           batch = basename(dirname(idat_red)))

# Read idats
raw_intensity_data <- read.metharray(samplesheet$idat_basename, force=T)
# Add sample IDs
colnames(raw_intensity_data) <- samplesheet$sample

# Add annotation
annotation(raw_intensity_data) <- c(
  array = "IlluminaHumanMethylationEPIC",
  annotation = "ilm10b5.hg38"
)

# Read lists of probes to filter
Zhou_probes <- read.delim(Zhou_input,col.names = 'Probe', header = F)
CrossReactive_probes <- read.delim(CrossReactive_input)
Problematic_probes <- read.delim(Problematic_input, col.names = 'Probe')
Filter_probes <- unique(c(Zhou_probes$Probe, CrossReactive_probes$Probe, Problematic_probes$Probe))
#-------------------------------------------------------------------------------
# 2.1 Quality Control: filter out samples with >5% failed probes
#-------------------------------------------------------------------------------
# Calculate detection P values
detection_pvalues <- raw_intensity_data %>% detectionP()

# Identify samples with more than 5% failed probes (pvalue cutoff 0.01)
failed_samples <- colMeans(detection_pvalues > 0.01) > 0.05
if(any(failed_samples)){
    message("Failed samples: ", paste(colnames(raw_intensity_data)[failed_samples], collapse=", "))
    # Filter out failed samples
    raw_intensity_data <- raw_intensity_data[, !failed_samples]
    detection_pvalues <- detection_pvalues[,!failed_samples]
}


#-------------------------------------------------------------------------------
# 2.2 Normalization: Perform Noob normalization
#-------------------------------------------------------------------------------
normalized_data <- preprocessNoob(raw_intensity_data)

#-------------------------------------------------------------------------------
# 2.3 Filter probes
#-------------------------------------------------------------------------------
# Keep only probes that succeeded in all samples
keep_probes <- rowSums(detection_pvalues < 0.01) == ncol(normalized_data)
normalized_data <- normalized_data[keep_probes, ]

# Filter problematic probes 
normalized_data <- normalized_data[
  !rownames(normalized_data) %in% Filter_probes,
  ]

# Add genomic coordinates
normalized_data <- mapToGenome(normalized_data)

# Remove SNP probes
normalized_data <- dropLociWithSnps(normalized_data)

# remove sex chromosomes
annotation_df <- getAnnotation(normalized_data)
normalized_data <- normalized_data[!(annotation_df$chr %in% c("chrX","chrY")), ]

message("Remaining probes after filtering: ", nrow(normalized_data))

#-------------------------------------------------------------------------------
# 2.4 Extract data
#-------------------------------------------------------------------------------
# Fetch methylation data
beta_values <- getBeta(normalized_data)
m_values <- getM(normalized_data)

# Fetch sample metadata
sample_metadata <- samplesheet %>% select(patient,sample,batch)

# fetch probe metadata
probe_metadata <- as.data.frame(getAnnotation(normalized_data)) %>%
    select(-c(ProbeSeqA,Forward_Sequence,SourceSeq)) %>%
    # Replace NA with empty string from compatibility
    mutate(across(everything(), ~ ifelse(is.na(.), "", .)))

#-------------------------------------------------------------------------------
# 3.0 Create AnnData object
#-------------------------------------------------------------------------------
adata <- anndata::AnnData(
  X = t(m_values),         
  obs = sample_metadata,        
  var = probe_metadata,
)
adata$layers[['beta']] <- t(beta_values)
#-------------------------------------------------------------------------------
# 4.0 Write anndata to .h5ad 
#-------------------------------------------------------------------------------
write_h5ad(adata, output)

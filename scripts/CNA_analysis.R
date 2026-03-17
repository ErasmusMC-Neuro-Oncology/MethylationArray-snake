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
    output <-  snakemake@output[[1]]
}else{
    input <- '/home/jurriaan/Projects/MINT/data/samplesheets/samplesheet_methylation.csv'
    input_reference <- '/home/jurriaan/Projects/Capper_Methylation/data/samplesheet.csv'
    output <- ''
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
array <- strsplit(annotation(raw_intensity_data)[1],'IlluminaHumanMethylation')[[1]][2]
#-------------------------------------------------------------------------------
# 3.1 CNV analysis
#-------------------------------------------------------------------------------
# Fetch annotations and subset probes
anno <- CNV.create_anno(array_type = array)
anno@probes <- anno@probes[names(anno@probes) %in% names(minfi::getLocations(IlluminaHumanMethylationEPICanno.ilm10b4.hg19::IlluminaHumanMethylationEPICanno.ilm10b4.hg19))]

# create CNV object
CNV_object <- CNV.load(normalized_data)


# Estimate CNVs and save objects in list
tumor_samples <- which(samplesheet$group == 'query')
cnv_list <- lapply(tumor_samples, function(i) {
    CNV.fit(CNV_object[i,],CNV_object[samplesheet$group == 'reference',],anno=anno )
})


CNV_MINT20 <- CNV.fit(CNV_object['MINT20_tumor1',],CNV_object[samplesheet$group == 'reference',],anno=anno )
binned <- CNV.bin(CNV_MINT20)
segmented <- CNV.segment(binned)

pdf('CNA_profile_MINT20_EPIC.pdf', width = 6 , height = 5)
CNV.genomeplot(segmented)
dev.off()
samplesheet %>% filter(sample == 'MINT20_tumor1')
# Perform binning and segmentation
cnv_list <- lapply(cnv_list, CNV.bin)
cnv_list <- lapply(cnv_list, CNV.segment)

plot()
CNV.fit(CNV_object[1,],CNV_object[samplesheet$group == 'reference',],anno=anno )


normalized_data


# Estimate CNVs per sample and store results in list
tumor_samples <- which(samplesheet$group == 'query')
cnv_list <- lapply(tumor_samples, function(i) {
  CNV.fit(
    CNV_object[, i],
    reference = CNV_object[samplesheet$group == 'reference',],
    anno = anno
  )
})

# Perform binning
cnv_list <- lapply(cnv_list, CNV.bin)

# Segment CNVs
cnv_list <- lapply(cnv_list, CNV.segment)







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
# 2.5 Extract data
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

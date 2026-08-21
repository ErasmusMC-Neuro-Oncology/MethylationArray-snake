#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Calculate_CGC.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Calculate Continuous Grading Coefficient (CGC)
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
#  04-06-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(minfi))
suppressMessages(library(dplyr))
library(IlluminaHumanMethylationEPICv2manifest)
library(IlluminaHumanMethylationEPICv2anno.20a1.hg38)

# Use correct reticulate environment
reticulate::use_condaenv(Sys.getenv("CONDA_PREFIX"), required = TRUE)
#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input<- snakemake@input[[1]]
    classifier <- snakemake@params[['classifier']]
    utils <- snakemake@params[['utils']]
    output <- snakemake@output[[1]]
}else{
    input <- '/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/samplesheets/Samplesheet_Methylation.csv'
    utils <- '/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/sslowgrade/workflows/MethylationArray-snake/scripts/anndata_utils.R'
    classifier <- '/home/jurriaan/mnt/BIGR_home/Resources/Continuous_Grading_Classifier/assets/CGC-Psi_predictor_probe_based_lm_v1.0_epic.Rds'
    output <- '/home/jurriaan/mnt/BIGR_home/SSLOWGRADE/output/Methylation/results/CGC.txt'
}
#-------------------------------------------------------------------------------
# 1.1 Read data
#-------------------------------------------------------------------------------
# Read .idat files
samplesheet <- read.delim(input , sep = ',')  %>%
    mutate(idat_basename = gsub("_Red.idat$", "", idat_red),
           batch = basename(dirname(idat_red)))

RGSet <- read.metharray(samplesheet$idat_basename, force=T, verbose = T)

if(RGSet@annotation[[1]] == 'Unknown'){
    annotation(RGSet) <- c(array = "IlluminaHumanMethylationEPICv2", 
                        annotation = "20a1.hg38")
}

#-------------------------------------------------------------------------------
# 2.1 Fetch M values
#-------------------------------------------------------------------------------
proc <- preprocessNoob(RGSet, offset = 0, dyeCorr = T, verbose = TRUE, dyeMethod="single")  #dyeMethod="reference"

mvalue <- ratioConvert(proc, what = "M") |>
  assays() |>
  purrr::pluck('listData') |>
  purrr::pluck("M") |>
  data.table::as.data.table(keep.rownames = "probe_id")
#-------------------------------------------------------------------------------
# 2.2 Predict CGC
#-------------------------------------------------------------------------------
classifier <- readRDS(classifier)

# acquire the exact same m-values
data <- mvalue |> 
  tibble::column_to_rownames('probe_id') |> 
  t() |> 
  as.data.frame() |> 
  dplyr::select(rownames(predictor$beta)) |> 
  as.matrix()


# apply lm to the data
out <- glmnet::predict.glmnet(predictor, data) |> 
  as.data.frame() |> 
  dplyr::rename(`CGCψ` = 1) |> 
  dplyr::rename_with(.fn = ~ paste0(., suffix), .cols = c('CGCψ'))

#-------------------------------------------------------------------------------
# 3.1 Join an write to file
#-------------------------------------------------------------------------------
write.table(out,output, sep = '\t',quote = F, row.names = F)

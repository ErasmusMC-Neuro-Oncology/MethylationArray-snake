#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot_SampleData.R
#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
# Create SampleData plots
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
#  21-08-2026: File creation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 0.1  Load packages
#-------------------------------------------------------------------------------
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))

#-------------------------------------------------------------------------------
# 0.2 Parse command line arguments
#-------------------------------------------------------------------------------
if(exists("snakemake")){
    input_SampleData <- snakemake@input[["SampleData"]]
    output <- snakemake@output[["CGC"]]
}else{
    input_SampleData <- "~/mnt/BIGR_home/SSLOWGRADE/output/Methylation/sampledata/SampleData_Methylation.txt"
    output <- "~/mnt/BIGR_home/SSLOWGRADE/output/Methylation/plots/CGC_primary_recurrence.pdf"
}
#-------------------------------------------------------------------------------
# 1.1 Read data 
#-------------------------------------------------------------------------------
# Read datasets
SampleData <- read.delim(input_SampleData)

#-------------------------------------------------------------------------------
# 1.2 Plot data
#-------------------------------------------------------------------------------
SampleData <- SampleData %>%
    mutate(patient = gsub('_R1_tumor1|_tumor1','',sample),
           Surgery = factor(ifelse(grepl('R1',sample),'primary','recurrence'), levels = c("primary", "recurrence")),
           x_jit   = as.numeric(Surgery) + runif(n(), -0.12, 0.12))


# Plot paired boxplot CGC
pdf(output, height = 4, width = 4)
SampleData %>% filter(!is.na(CGC_curated)) %>%
ggplot(aes(x = x_jit, y = CGC_curated)) +
  geom_boxplot(aes(x = as.numeric(Surgery), group = Surgery),
               outlier.shape = NA, width = 0.5, fill = NA) +
  geom_line(aes(group = patient), colour = "grey50",
            linewidth = 0.3, alpha = 0.6) +
  geom_point(aes(color = predicted_class), size = 1.5) +
  scale_x_continuous(breaks = 1:2, labels = levels(df$Surgery),
                     expand = expansion(add = 0.6)) +
  labs(x = "", y = "CGC (curated)",color = '') +
  theme_classic() + theme(legend.position = "top") 
dev.off()

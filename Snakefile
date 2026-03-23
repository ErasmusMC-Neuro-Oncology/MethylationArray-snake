configfile: "config.yaml"
from datetime import datetime
#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare wildcards and variables
output_dir = config["all"]["output_dir"]
datasets = config['all']['datasets']
#-------------------------------------------------------------------------------------------------------------------
# 0.2 specify target rules
rule all:
    input:
        expand(output_dir + "methylation/{dataset}/methylation_data.h5ad",dataset = datasets)

#+++++++++++++++++++++++++++++++++++++++++ 1. PREPROCESS IDAT FILES  +++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Perform QC, normalization and compute Beta/M-values
rule Preprocess_idat:
    input:
        lambda wildcards: config['samplesheet'][wildcards.dataset]
    output:
        Mset = output_dir + "methylation/{dataset}/methylation_object.Rds",
        adata = output_dir + "methylation/{dataset}/methylation_data.h5ad"
    conda:
        'envs/minfi.yaml'
    params:
        Zhou_probes = config['EPIC']['Zhou'],
        CrossReactive_probes = config['EPIC']['CrossReactive'],
        Problematic_probes = config['EPIC']['Problematic']
    threads: 2
    resources:
        mem_mb=10000
    script:
        "scripts/Preprocess_idat.R"

#-------------------------------------------------------------------------------------------------------------------
# 1.2 Perform CNA analysis, use Pai et al normals as a reference
rule CNA_analysis:
    input:
        query = output_dir + "methylation/MINT/methylation_object.Rds",
        reference =  output_dir + "methylation/Pai/methylation_object.Rds",
    output:
        Segmented = output_dir + "CNAs/MINT/Segmented_CNAs_MINT.txt",
        Profile_dir = directory(output_dir + 'CNAs/MINT/plots/')
    conda:
        'envs/minfi.yaml'
    threads: 2
    resources:
        mem_mb=10000
    script:
        "scripts/CNA_analysis.R"


#+++++++++++++++++++++++++++++++++++++++++ 2. ESTIMATE TUMOR PURTIY  +++++++++++++++++++++++++++++++++++++++++++++
# 2.1 Estimate tumor purtiy using RF_purity and InfiniumPurify
rule Estimate_tumor_purity:
    input:
        output_dir + "methylation/methylation_data_{dataset}.h5ad",
    output:
        output_dir + "results/Tumor_purities.txt"
    conda:
        'envs/minfi.yaml'
    threads: 2
    resources:
        mem_mb=10000
    script:
        "scripts/Estimate_tumor_purtiy.R"


        

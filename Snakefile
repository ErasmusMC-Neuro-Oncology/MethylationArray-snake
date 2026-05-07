configfile: "config.yaml"
import glob as _glob
import os
from datetime import datetime
#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare wildcards and variables
output_dir = config["all"]["output_dir"]
datasets = config['all']['datasets']
normalizations = config['all']['normalizations']
reference_datasets = config['all']['reference_datasets']
cna_reference = config['all']['cna_reference']
query_datasets = [d for d in datasets if d not in reference_datasets]
download_datasets = config['all'].get('download_datasets', [])

wildcard_constraints:
    dataset = "[^/]+"

def get_samplesheet(wildcards):
    if wildcards.dataset in download_datasets:
        checkpoints.download_idat.get(dataset=wildcards.dataset)
        return output_dir + "idat/{}/samplesheet.csv".format(wildcards.dataset)
    return config['samplesheet'][wildcards.dataset]

def get_idat_inputs(wildcards):
    if wildcards.dataset in download_datasets:
        checkpoints.download_idat.get(dataset=wildcards.dataset)
        idat_dir = output_dir + "idat/{}/".format(wildcards.dataset)
        return sorted(_glob.glob(os.path.join(idat_dir, "*.idat")))
    return []

#-------------------------------------------------------------------------------------------------------------------
# 0.2 specify target rules
rule all:
    input:
        expand(output_dir + "methylation/{dataset}/methylation_data_{norm}.h5ad",dataset = datasets, norm = normalizations),
        expand(output_dir + "CNAs/{dataset}/Segmented_CNAs_{dataset}.txt", dataset = query_datasets),
        output_dir + "results/Tumor_purities.txt"


#+++++++++++++++++++++++++++++++++++++++++ 1. PREPROCESS IDAT FILES  +++++++++++++++++++++++++++++++++++++++++++++
# 1.0 Download IDAT files from URLs defined in samplesheet
checkpoint download_idat:
    input:
        lambda wildcards: config['samplesheet'][wildcards.dataset]
    output:
        directory(output_dir + "idat/{dataset}/"),
        output_dir + "idat/{dataset}/.download_complete"
    script:
        "scripts/download_idat.py"

#-------------------------------------------------------------------------------------------------------------------
# 1.1 Perform QC, normalization and compute Beta/M-values
rule Preprocess_idat:
    input:
        samplesheet = get_samplesheet,
        idats       = get_idat_inputs
    output:
        Mset = output_dir + "methylation/{dataset}/methylation_data_{norm}.Rds",
        adata = output_dir + "methylation/{dataset}/methylation_data_{norm}.h5ad"
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
# 1.2 Perform CNA analysis; reference dataset fixed via config['all']['cna_reference']
rule CNA_analysis:
    input:
        query = output_dir + "methylation/{dataset}/methylation_data_noob.Rds",
        reference = output_dir + "methylation/" + cna_reference + "/methylation_data_noob.Rds"
    output:
        Segmented = output_dir + "CNAs/{dataset}/Segmented_CNAs_{dataset}.txt",
        Profile_dir = directory(output_dir + "CNAs/{dataset}/plots/")
    conda:
        'envs/minfi.yaml'
    threads: 2
    resources:
        mem_mb=10000
    log:
        output_dir + "CNAs/{dataset}/CNA_analysis.log"
    script:
        "scripts/CNA_analysis.R"


#+++++++++++++++++++++++++++++++++++++++++ 2. ESTIMATE TUMOR PURTIY  +++++++++++++++++++++++++++++++++++++++++++++
# 2.1 Estimate tumor purtiy using RF_purity and InfiniumPurify
rule Estimate_tumor_purity:
    input:
        output_dir + "methylation/MINT/methylation_data_noob.h5ad"
    output:
        output_dir + "results/Tumor_purities.txt"
    conda:
        'envs/minfi.yaml'
    threads: 2
    resources:
        mem_mb=10000
    script:
        "scripts/Estimate_tumor_purity.R"



#+++++++++++++++++++++++++++++++++++++++++ 3. CLASSIFICATION  +++++++++++++++++++++++++++++++++++++++++++++
# 3.1 Classify
rule Estimate_tumor_subtype:
    input:
        output_dir + "methylation/MINT/methylation_data_noob.h5ad"
    output:
        output_dir + "results/Tumor_purities.txt"
    conda:
        'envs/minfi.yaml'
    threads: 2
    resources:
        mem_mb=10000
    script:
        "scripts/Estimate_tumor_purity.R"

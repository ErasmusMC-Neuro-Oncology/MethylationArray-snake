configfile: "config.yaml"
from datetime import datetime
#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare wildcards and variables
output_dir = config["all"]["output_dir"]

#-------------------------------------------------------------------------------------------------------------------
# 0.2 specify target rules
rule all:
    input:
        output_dir + 'sampledata/SampleData_Methylation.txt'

#+++++++++++++++++++++++++++++++++++++++++ 1. PREPROCESS IDAT FILES  +++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Perform QC, normalization and compute Beta/M-values
rule Preprocess_idat:
    input:
        config['all']['samplesheet']
    output:
        Mset = output_dir + "methylation_object.Rds",
        adata = output_dir + "methylation_data.h5ad"
    conda:
        'envs/minfi.yaml'
    params:
        Zhou_probes = config['EPIC']['Zhou'],
        CrossReactive_probes = config['EPIC']['CrossReactive'],
        Problematic_probes = config['EPIC']['Problematic']
    threads: 2
    resources:
        mem_mb=100000,
        gpu = 0
    script:
        "scripts/Preprocess_idat.R"

#-------------------------------------------------------------------------------------------------------------------
# 1.2 Perform CNA analysis, use Pai et al normals as a reference
rule CNA_analysis:
    input:
        query = output_dir + "methylation_object.Rds",
    output:
        Segmented = output_dir + "CNAs/Segmented_CNAs.txt",
        Profile_dir = directory(output_dir + 'CNAs/plots/')
    params:
        reference =  config['all']['reference'],
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
        output_dir + "methylation_data.h5ad"
    output:
        output_dir + "results/Tumor_purities.txt"
    params:
        utils = config['all']['utils']
    conda:
        'envs/minfi.yaml'
    threads: 2
    resources:
        mem_mb=10000
    script:
        "scripts/Estimate_tumor_purity.R"



#+++++++++++++++++++++++++++++++++++++++++ 3. CLASSIFICATION  +++++++++++++++++++++++++++++++++++++++++++++
# 3.1 Classify using pretrained model
rule Classify_samples:
    input:
        config['all']['samplesheet']
    output:
        output_dir + "results/Methylation_Classes.txt"
    params:
        classifier = config['classify']['classifier'],
        ba_coef = config['classify']['ba_coef'],
        material = config['classify']['material'],
        Rpreprocess = config['classify']['preprocess_script'],
        filter_dir = config['classify']['filter_dir'],
        CNA_data = config['classify']['CNA_data']
    conda:
        'envs/classify.yaml'
    resources:
        mem_mb=10000
    script:
        "scripts/Classify_samples.R"

        
#+++++++++++++++++++++++++++++++++++++++++ 4. CREATE SAMPLE DATA  +++++++++++++++++++++++++++++++++++++++++++++
# 3.1 Classify using pretrained model

rule Create_SampleData:
    input:
        Classes = output_dir + "results/Methylation_Classes.txt",
        purities = output_dir + "results/Tumor_purities.txt"
    output:
        output_dir + 'sampledata/SampleData_Methylation.txt'
    conda:
        "envs/R.yaml"
    script:
        'scripts/Create_SampleData.R'


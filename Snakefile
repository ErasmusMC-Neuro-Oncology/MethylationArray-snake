configfile: "config.yaml"
from datetime import datetime
#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare wildcards and variables
data_dir = config["all"]["data_dir"]
output_dir = config["all"]["output_dir"]
#-------------------------------------------------------------------------------------------------------------------
# 0.2 specify target rules
rule all:
    input:
        output_dir + "methylation/methylation_data.h5ad"

#+++++++++++++++++++++++++++++++++++++++++ 1. PREPROCESS IDAT FILES  +++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Perform QC, normalization and compute Beta/M-values
rule Preprocess_idat:
    input:
        config['all']['samplesheet']
    output:
        output_dir + "methylation/methylation_data.h5ad"
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
# 1.2 

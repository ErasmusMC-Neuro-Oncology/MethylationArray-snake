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
        "results/MethylArray/.done"

#+++++++++++++++++++++++++++++++++++++++++ 1. PREPROCESS IDAT FILES  +++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Calculate M-value
rule Preprocess_idat:
    input:
        config['all']['samplesheet']
    output:
        output_dir + "methylation/M_values.Rds"
    conda:
        'envs/minfi.yaml'
    threads: 2
    resources:
        mem_mb=10000
    script:
        "scripts/Preprocess_idat.R"

#-------------------------------------------------------------------------------------------------------------------
# 1.2

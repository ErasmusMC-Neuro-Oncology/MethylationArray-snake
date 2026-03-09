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

#+++++++++++++++++++++++++++++++++++++++++ 1 RUN METHYLARRAY  +++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Run nf-core MethylArray
rule MethylArray:
    input:
        "../MINT/data/samplesheets/samplesheet_methylation.csv"
    output:
        "results/MethylArray/.done"
    threads: 2
    resources:
        mem_mb=1000
    conda:
        "envs/nextflow.yaml"
    log:
        "logs/MethylArray/nextflow_"+datetime.now().strftime("%Y_%m_%d_%H:%M:%S")+".log"
    params:
        genome="hg38",
        profile="singularity",
        outdir = output_dir + 'MethylArray'
    shell:
        """
        nextflow -log {log} run  nf-core/methylarray -r eb5fb7d \
            -profile {params.profile} \
            --input {input} \
            --outdir {params.outdir} \
            --bs_genome_version {params.genome} \
            --max_cpus {threads} \
            --max_memory '{resources.mem_mb} MB'
            -resume 
            
        touch {output}
        """

#-------------------------------------------------------------------------------------------------------------------
# 1.2

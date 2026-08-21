configfile: "config.yaml"
from datetime import datetime
import pandas as pd
#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare wildcards and variables
output_dir = config["all"]["output_dir"]
Samples = pd.read_csv(config['all']['samplesheet'])['patient']
#-------------------------------------------------------------------------------------------------------------------
# 0.2 specify target rules
rule all:
    input:
        output_dir + 'sampledata/SampleData_Methylation.txt',
        output_dir + 'CGC/CGC_psi.tsv'

#+++++++++++++++++++++++++++++++++++++++++ 1. PREPROCESS IDAT FILES  +++++++++++++++++++++++++++++++++++++++++++++
# 1.1 Perform QC, normalization, compute Beta/M-values and compute segments
rule Preprocess_idat:
    input:
        config['all']['samplesheet']
    output:
        adata = output_dir + "methylation_data.h5ad"
    conda:
        'envs/methylation.yaml'
    params:
        Zhou_probes = config['EPIC']['Zhou'],
        CrossReactive_probes = config['EPIC']['CrossReactive'],
        Problematic_probes = config['EPIC']['Problematic']
    threads: 20
    resources:
        mem='100GB',
        gpu = 0
    script:
        "scripts/Preprocess_idat.R"

#+++++++++++++++++++++++++++++++++++++++++ 2. CLASSIFY SAMPLES  +++++++++++++++++++++++++++++++++++++++++++++
# 2.1 Classify samples using pretrained classifier
rule crossNN:
    input:
        adata = output_dir + "methylation_data.h5ad"
    output:
        predictions = output_dir + 'crossNN/crossNN_predictions.tsv',
        scores      = output_dir + 'crossNN/crossNN_scores.tsv',
        coverage    = output_dir + 'crossNN/crossNN_feature_coverage.tsv'
    params:
        snake_dir    = config['all']['snake_dir'],
        model_dir    = config['crossNN']['model_dir'],
        repo_url     = config['crossNN']['repo_url'],
        commit       = config['crossNN']['commit'],
        weights      = config['crossNN']['weights'],
        score_cutoff = config['crossNN']['score_cutoff'],
        min_features = config['crossNN']['min_features']
    conda:
        "envs/crossNN.yaml"
    threads: 4
    resources:
        mem_mb=32000,
        runtime='1h',
        gpu=0
    shell:
        """
        # Fetch the crossNN model (skipped if already present)
        if [ ! -d {params.model_dir}/.git ]; then
            rm -rf {params.model_dir}
            git clone {params.repo_url} {params.model_dir}
            cd {params.model_dir}
            git checkout {params.commit}
            git lfs install --local && git lfs pull
            cd -
        fi

        # Run crossNN
        python3 {params.snake_dir}/scripts/Run_crossNN.py \
        -i {input.adata} \
        -m {params.model_dir} \
        -w {params.weights} \
        -s {params.score_cutoff} \
        -n {params.min_features} \
        -t {threads} \
        -o_pred {output.predictions} \
        -o_scores {output.scores} \
        -o_coverage {output.coverage}
        """

# 2.2 Compute CGC-Psi from the harmonized methylation data
rule CGC:
    input:
        adata = output_dir + "methylation_data.h5ad"
    output:
        cgc = output_dir + 'CGC/CGC_psi.tsv'
    params:
        predictor     = config['CGC']['predictor'],
        predictor_url = config['CGC']['predictor_url'],
        layer         = config['CGC']['layer'],
        value_type    = config['CGC']['value_type']
    conda:
        "envs/CGC.yaml"
    threads: 1
    resources:
        mem_mb=32000,
        runtime='1h',
        gpu=0
    script:
        "scripts/Run_CGC.R"


#+++++++++++++++++++++++++++++++++++++++ 3. COPY NUMBER PLOTS +++++++++++++++++++++++++++++++++++++++
# 3.1 Draw genome-wide segmented CNV plots, one per sample
rule CNV_plots:
    input:
        adata = output_dir + "methylation_data.h5ad"
    output:
        plots = expand(output_dir + 'CNV/{sample}_CNV.pdf', sample = Samples)
    params:
        outdir  = output_dir + 'CNV',
        samples = list(Samples),
        suffix  = '_CNV.pdf',
        ymin    = config['CNV']['ymin'],
        ymax    = config['CNV']['ymax'],
        width   = config['CNV']['width'],
        height  = config['CNV']['height']
    conda:
        "envs/CNV.yaml"
    threads: 1
    resources:
        mem_mb=64000,
        runtime='2h',
        gpu=0
    script:
        "scripts/Plot_CNV.R"

#+++++++++++++++++++++++++++++++++++++++ 4. CREATE SAMPLEDATA +++++++++++++++++++++++++++++++++++++++
# 4.1 Draw genome-wide segmented CNV plots, one per sample
rule Create_SampleData:
    input:
        predictions = output_dir + 'crossNN/crossNN_predictions.tsv',
        scores = output_dir + 'crossNN/crossNN_scores.tsv',
        cgc = output_dir + 'CGC/CGC_psi.tsv'
    output:
        output_dir + 'sampledata/SampleData_Methylation.txt'
    conda:
        "envs/R.yaml"
    script:
        'scripts/Create_SampleData.R'

#++++++++++++++++++++++++++++++++++++++++++++++++ 5 PLOT SAMPLE DATA +++++++++++++++++++++++++++++++++++++++++++++++++++++
# Plot SampleData
rule Plot_SampleData:
    input:
        SampleData = output_dir + 'sampledata/SampleData_Methylation.txt'
    output:
        CGC = output_dir + 'plots/CGC_primary_recurrence.pdf'
    conda:
        "envs/R.yaml"
    script:
        'scripts/Plot_SampleData.R'
        


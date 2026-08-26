configfile: "config.yaml"
from datetime import datetime
#+++++++++++++++++++++++++++++++++++++++ 0 PREPARE WILDCARDS AND TARGET ++++++++++++++++++++++++++++++++++++++++++++
# 0.1 Prepare wildcards and variables
output_dir = config["all"]["output_dir"]
data_dir =  config["all"]["data_dir"]
#-------------------------------------------------------------------------------------------------------------------
# 0.2 specify target rules
rule all:
    input:
        output_dir + 'sampledata/SampleData_Methylation.txt'

#+++++++++++++++++++++++++++++++++ 1. HARMONIZE + BATCH CORRECT ++++++++++++++++++++++++++++++++++
rule Preprocess_methylation:
    input:
        samplesheet_mint = config['all']['samplesheet']
    output:
        adata = output_dir + "harmonized/methylation_harmonized.h5ad"
    params:
        out_dir = output_dir + 'harmonized',
        idat_dir_Capper = config['harmonize']['idat_dir_Capper'],
        idat_dir_Lucas = config['harmonize']['idat_dir_Lucas'],
        samplesheet_Capper = config['harmonize']['samplesheet_Capper'],
        samplesheet_Sturm = config['harmonize']['samplesheet_Sturm'],
        samplesheet_Lucas = config['harmonize']['samplesheet_Lucas'],
        classes_csv = config['harmonize']['classes_csv'],
        Zhou_probes = config['EPIC']['Zhou'],
        CrossReactive_probes = config['EPIC']['CrossReactive'],
        Problematic_probes = config['EPIC']['Problematic'],
    conda:
        "envs/methylation.yaml"
    threads: 20
    resources:
        mem_mb=400000,
        runtime='24h',
        gpu=0
    script:
        "scripts/Preprocess_methylation.R"

#++++++++++++++++++++++++++++++++++++ 2. PCA / t-SNE EMBEDDINGS ++++++++++++++++++++++++++++++++++
rule Embeddings_methylation:
    input:
        adata = output_dir + "harmonized/methylation_harmonized.h5ad"
    output:
        embeddings_all  = output_dir + 'embeddings/Embeddings_and_subtypes.csv',
        embeddings_selected = output_dir + 'embeddings/Embeddings_and_subtypes_selected.csv',
        tsne_family= output_dir + 'embeddings/tSNE_AllGliomas_w_batch_correction.pdf',
        tsne_zoom = output_dir + 'embeddings/tSNE_AllGliomas_zoom_tsne_-20_0.pdf',
        tsne_pathology = output_dir + 'embeddings/tSNE_AllGliomas_Histological_subtype.pdf',
        pca_before = output_dir + 'embeddings/PCA_before_ComBat.png',
        pca_after = output_dir + 'embeddings/PCA_after_ComBat.png'
    params:
        out_dir           = output_dir + 'harmonized',
    conda:
        "envs/methylation_embeddings.yaml"
    threads: 16
    resources:
        mem_mb=200000,
        runtime='8h',
        gpu=0
    script:
        "scripts/Embeddings_Methylation.R"
        
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


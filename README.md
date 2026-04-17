# MethylationArray-snake

## Overview

This pipeline performs analysis of **Illumina DNA methylation array data** using a Snakemake workflow. It is designed for processing raw IDAT files from EPIC / Infinium methylation arrays.

The workflow combines:

* **minfi** for IDAT import, quality control, preprocessing and normalization
* Probe filtering for problematic, cross-reactive and SNP-associated probes
* Conversion to `.h5ad` format for downstream analysis in Python / Scanpy workflows
* Copy number analysis derived from methylation intensities
* Tumor purity estimation using methylation-based models

The main goal is to generate high-quality normalized methylation data together with copy number profiles and tumor purity estimates.

---

## Workflow Summary

### 1. Sample Preparation

Input samplesheets are defined in `config.yaml` for each dataset. This script has to be edited by the user.
 
### 2. Preprocessing of IDAT Files

For each configured dataset:

* Raw IDAT files are imported
* Quality control is performed
* Signal intensities are normalized
* Problematic probes are removed
* Beta-values / M-values are computed
* Data are exported as R objects and `.h5ad` files

### 3. Copy Number Analysis

Using methylation intensity data:

* CNA profiles are inferred from array signal intensities
* The `MINT` cohort is analyzed using the `Pai` cohort as reference normals
* Segmented CNA profiles and plots are generated

### 4. Tumor Purity Estimation

Tumor purity is estimated using:

* **RF_purity**
* **InfiniumPurify**

### 5. Reporting / Export

Processed methylation objects, CNA profiles and purity tables are written to disk.

---

## Software Requirements

* **Snakemake 9.17.3**
* Conda / Mamba

---

## Running the Pipeline

```bash
snakemake --cores 1 --use-conda
```

---

## Main Output Files

### Normalized Methylation Data

```text
methylation/{dataset}/methylation_data_{norm}.h5ad
```

Primary output containing processed methylation data for downstream analysis.

### R Methylation Objects

```text
methylation/{dataset}/methylation_object_{norm}.Rds
```

Contains normalized minfi objects for R-based downstream analyses.

### Copy Number Output

```text
CNAs/MINT/Segmented_CNAs_MINT.txt
```

Primary CNA segmentation output derived from methylation array signal intensities.

### Tumor Purity Estimates

```text
results/Tumor_purities.txt
```

Per-sample tumor purity estimates generated from methylation data using RFpurity and InfiniumPurify.

### Additional Outputs

#### CNA Plots

```text
CNAs/MINT/plots/
```

Contains per-sample CNA profile visualizations.

---

## Directory Structure

```text
output/
├── methylation/
├── CNAs/
└── results/
```

---

## Notes on Illumina Array Processing

This pipeline supports Illumina methylation array workflows where preprocessing quality is essential.

By default, problematic probes can be removed using supplied reference lists:

* Zhou probes
* Cross-reactive probes
* Other problematic probes

Normalization methods are defined in `config.yaml`.

---

#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "${SCRIPT_DIR}")"

cd "${PROJECT_DIR}"

source ~/.zshrc
conda activate snakemake

snakemake \
    --use-conda \
    --cores 2 \
    'output/CNAs/test/Segmented_CNAs_test.txt'

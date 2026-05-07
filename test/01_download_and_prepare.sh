#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
IDAT_DIR="${SCRIPT_DIR}/idat"

# Download IDAT files
echo "Downloading IDAT files..."
wget -q -P "${IDAT_DIR}" \
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8997664&format=file&file=GSM8997664%5F206467010068%5FR06C01%5FGrn%2Eidat%2Egz" \
    -O "${IDAT_DIR}/GSM8997664_206467010068_R06C01_Grn.idat.gz"

wget -q -P "${IDAT_DIR}" \
    "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSM8997664&format=file&file=GSM8997664%5F206467010068%5FR06C01%5FRed%2Eidat%2Egz" \
    -O "${IDAT_DIR}/GSM8997664_206467010068_R06C01_Red.idat.gz"

# Unzip
echo "Extracting..."
gunzip -f "${IDAT_DIR}"/*.idat.gz

# Create samplesheet
SAMPLESHEET="${SCRIPT_DIR}/samplesheet.csv"
echo "sample,idat_red" > "${SAMPLESHEET}"
echo "GSM8997664,${IDAT_DIR}/GSM8997664_206467010068_R06C01_Red.idat" >> "${SAMPLESHEET}"

echo "Done. Samplesheet written to ${SAMPLESHEET}"

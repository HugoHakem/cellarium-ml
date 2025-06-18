#!/usr/bin/env bash
set -euo pipefail

# Define paths
DATA_DIR="fine_tune/datasets/genome"
ZIP_PATH="${DATA_DIR}/human_genes.zip"
TSV_PATH="${DATA_DIR}/human_gene_table.tsv"
EXTRACT_DIR="${DATA_DIR}/ncbi_dataset"
JSONL_PATH="${EXTRACT_DIR}/data/data_report.jsonl"

# Create target directory
mkdir -p "${DATA_DIR}"

echo "Downloading human gene metadata to ${ZIP_PATH}..."
datasets download gene taxon 9606 --filename "${ZIP_PATH}"

echo "Unzipping dataset..."
unzip -o "${ZIP_PATH}" -d "${DATA_DIR}"

echo "Extracting gene table to ${TSV_PATH}..."
dataformat tsv gene \
  --fields symbol,gene-type,name-authority,ensembl-geneids,synonyms \
  < "${JSONL_PATH}" \
  > "${TSV_PATH}"

# echo "Cleaning up..."
# rm -f "${ZIP_PATH}"
# rm -rf "${EXTRACT_DIR}"

echo "✅ Done. Gene table available at: ${TSV_PATH}"

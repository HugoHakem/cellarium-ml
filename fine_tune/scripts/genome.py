import gzip
import re
from pathlib import Path
from urllib import request

import pandas as pd
from jsonargparse import CLI

# Adaptation of https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-ref-build-steps#ref-2020-a

def download_ref_genome(
    data_path: Path=Path("./fine_tune/datasets/genome"),
    comp_name: str="gencode.v32.annotation.gtf.gz",
    out_name: str="refdata-gex-GRCh38-2020-A.csv"
):
    """
    Downloads and processes a GENCODE GTF annotation file to extract gene IDs and gene names.

    If the processed CSV file does not exist locally, this function downloads the GTF file,
    extracts gene_id and gene_name pairs for all genes, and saves them as a CSV file.
    If the CSV already exists, the function does nothing.

    Args:
        data_path (str or Path, optional): Directory to store the downloaded and processed files.
        comp_name (str, optional): Name for the downloaded GTF file.
        out_name (str, optional): Name for the output CSV file.

    Returns:
        None
    """
    if isinstance(data_path, str):
        data_path = Path(data_path)

    # Set paths
    gtf_url = "http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_32/gencode.v32.primary_assembly.annotation.gtf.gz"

    # Download GTF if not already downloaded
    if not (out_path := data_path.joinpath(out_name)).exists():
        if not (comp_path := data_path.joinpath(comp_name)).exists():
            comp_path.parent.mkdir(parents=True, exist_ok=True)
            print("Downloading GTF...")
            request.urlretrieve(gtf_url, comp_path)
            print("Done.")

        # Extract gene_id and gene_name
        genes = []
        pattern = re.compile(r'gene_id "([^"]+)";.*?gene_name "([^"]+)";')

        with gzip.open(comp_path, "rt") as f:
            print(f"Extracting {out_path}")
            for line in f:
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                if fields[2] == "gene":
                    attr = fields[8]
                    match = pattern.search(attr)
                    if match:
                        gene_id, gene_name = match.groups()
                        genes.append((gene_id.split('.')[0], gene_name))  # Strip version

        # Save to CSV
        df = pd.DataFrame(genes, columns=["ensemble_id", "symbol"]).drop_duplicates()
        df.to_csv(out_path := data_path.joinpath(out_name), index=False)

        print(f"Saved {len(df)} gene entries to {out_path}")
    else:
        print(f"Found local copy of {out_path}")

if __name__ == "__main__":
    CLI(download_ref_genome)

import gzip
import re
from pathlib import Path
from typing import Optional, Union
from urllib import request

import pandas as pd
from jsonargparse import CLI
from tqdm import tqdm

# Adaptation of https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-ref-build-steps#ref-2020-a

def build_ref_genome(
    data_path: Path=Path("./fine_tune/datasets/genome"),
    out_name: Optional[str]=None,
    version: Union[int, str]=32
) -> None:
    """
    Downloads and processes a GENCODE GTF annotation file to extract gene IDs and gene names.

    If the processed CSV file does not exist locally, this function downloads the GTF file,
    extracts gene_id and gene_name pairs for all genes, and saves them as a CSV file.
    If the CSV already exists, the function does nothing.

    Args:
        data_path (str or Path, optional): Directory to store the downloaded and processed files.
        out_name (str, optional): Name for the output CSV file.
            If not provided, the same as the decompressed gtf file.
        version (int or str): The release version of the gencode annotation

    Returns:
        None
    """
    if isinstance(data_path, str):
        data_path = Path(data_path)

    # Set paths
    gtf_url = ("http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
               f"release_{version}/gencode.v{version}.primary_assembly.annotation.gtf.gz")
    tmp_name = gtf_url.split("/")[-1]
    out_name = out_name or ".".join(tmp_name.split(".")[:-2]) + ".csv"

    # Download GTF if not already downloaded
    if (out_path := data_path.joinpath(out_name)).exists():
        print(f"Found local copy of {out_path}")
    else:
        if not (tmp_path := data_path.joinpath(tmp_name)).exists():
            tmp_path.parent.mkdir(parents=True, exist_ok=True)
            print(f"Downloading Gencode GTF file v{version}...")
            request.urlretrieve(gtf_url, tmp_path)
            print("Done.")

        # Extract gene_id and gene_name and filter out genes
        biotypes = {
            "protein_coding", "lncRNA",
            "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", "IG_V_gene",
            "IG_V_pseudogene", "IG_J_pseudogene", "IG_C_pseudogene",
            "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene",
            "TR_V_pseudogene", "TR_J_pseudogene"
        }
        gene_allowlist = set()
        gene_info = {}

        with gzip.open(tmp_path, "rt") as gtf:
            print(f"Extracting {tmp_path}")
            for line in tqdm(gtf, desc="Processing GTF", unit=" lines"):
                if line.startswith("#"):
                    continue
                fields = line.strip().split("\t")
                feature_type = fields[2]
                attributes_field = fields[8]

                # Parse attributes
                attr_pairs = re.findall(r'(\S+) "([^"]+)"', attributes_field)
                attrs = dict(attr_pairs)

                if feature_type == "transcript":
                    if attrs.get("gene_type") not in biotypes:
                        continue
                    if attrs.get("transcript_type") not in biotypes:
                        continue
                    tags = [v for k, v in attr_pairs if k == "tag"]
                    if "readthrough_transcript" in tags or "PAR" in tags:
                        continue
                    gene_id = attrs["gene_id"].split(".")[0]
                    gene_allowlist.add(gene_id)

                elif feature_type == "gene":
                    gene_id = attrs["gene_id"].split(".")[0]
                    gene_name = attrs.get("gene_name", "")
                    if gene_id not in gene_info:
                        gene_info[gene_id] = gene_name

        # Retain only allowed genes
        print(f"Retain only allowed genes: {len(gene_allowlist)}")
        filtered_gene_info = {
            gene_id: gene_info.get(gene_id, "") for gene_id in sorted(gene_allowlist)
        }
        # Save to CSV
        df = pd.DataFrame(filtered_gene_info.items(), columns=["ensemble_id", "symbol"])
        df.to_csv(out_path := data_path.joinpath(out_name), index=False)
        print(f"Saved {len(df)} gene entries to {out_path}")
        tmp_path.unlink()

if __name__ == "__main__":
    CLI(build_ref_genome)

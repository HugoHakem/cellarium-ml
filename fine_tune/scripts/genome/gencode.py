import gzip
import re
from pathlib import Path
from typing import Optional, Union
from urllib import request

import pandas as pd
from jsonargparse import CLI
from tqdm import tqdm

# Adaptation of
# https://www.10xgenomics.com/support/software/cell-ranger/downloads/cr-ref-build-steps#ref-2020-a
# Do not make use of cellranger
# Improvement on the mapping of gene_name

def download_url(
  url: str,
  path: Path,
  print_flag: Optional[str]=None
):
    """
    Download a file from a URL to a local path.

    Args:
        url (str): The URL to download from.
        path (Path): The local file path to save the download.
        print_flag (str, optional): A label to print during download for user feedback.

    Returns:
        None
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    print(f"Downloading {print_flag or 'file'}...")
    request.urlretrieve(url, path)
    print("Done.")


def get_hgnc_symbol(
    gene_df: pd.DataFrame,
    data_path: Path=Path("./fine_tune/datasets/genome"),
    hgnc_date: str="2021-03-01",
) -> pd.DataFrame:
    """
    Map gene annotations to standardized HGNC gene symbols.

    Downloads the HGNC table if not present, then maps Ensembl gene IDs and HGNC IDs
    in the input DataFrame to their approved HGNC gene symbols. Priority is given to
    Ensembl gene IDs, then HGNC IDs, and finally falls back to the original gene name
    if no mapping is found. Adds or updates the 'symbol' column in the DataFrame.

    Args:
        gene_df (pd.DataFrame): Input DataFrame with columns 'ensembl_id', 'hgnc_id', and 'symbol'.
        data_path (Path, optional): Directory to store or find the HGNC table.
        hgnc_date: (str, optional): date of the HGNC achive. Must figure in:
            https://www.genenames.org/download/archive/monthly/tsv/

    Returns:
        pd.DataFrame: Copy of the DataFrame with an updated 'symbol' column containing HGNC symbols.
    """

    # hgnc table retrieved from here: https://www.genenames.org/download/archive/monthly/tsv/:
    #   - 2021-03-01, the oldest archive and closest to the gencode v32 of Sept 2019.
    # For up to date info, one may want to investigate: https://www.genenames.org/download/custom/
    # with the selected col:
    #     "HGNC ID", "Approved symbol",
    #     "Ensembl ID(supplied by Ensembl)", "Ensembl gene ID"

    hgnc_url = (
        "https://storage.googleapis.com/public-download-files/hgnc/archive/"
        f"archive/monthly/tsv/hgnc_complete_set_{hgnc_date}.txt"
    )
    if not (tmp_hgnc_path := data_path.joinpath("tmp_hgnc.txt")).exists():
        download_url(hgnc_url, tmp_hgnc_path, print_flag="HGNC table of 2021-03-01")
    hgnc_df = pd.read_csv(tmp_hgnc_path, sep="\t")

    # Map ensembl_id and hgnc_id to hgnc symbol. Priority given to the ensembl_id, then hgnc, then gene_name
    print("Map HGNC symbol")
    hgnc_id_to_symbol = (
        hgnc_df[["hgnc_id", "symbol"]]
        .dropna()
        .set_index("hgnc_id")["symbol"]
        .to_dict()
    )
    ensembl_id_to_symbol = (
        hgnc_df[["ensembl_gene_id", "symbol"]]
        .dropna()
        .set_index("ensembl_gene_id")["symbol"]
        .to_dict()
    )
    gene_df_copy = gene_df.copy()
    gene_df_copy["symbol"] = gene_df_copy[["ensembl_id", "hgnc_id", "symbol"]].apply(
        lambda x: (
            ensembl_id_to_symbol.get(x["ensembl_id"]) or
            (hgnc_id_to_symbol.get(x["hgnc_id"]) or x["symbol"])
        ),
        axis=1
    )

    tmp_hgnc_path.unlink()
    print(f"Number of genes renamed with HGNC symbol: {(gene_df_copy['symbol'] != gene_df['symbol']).sum()}")
    return gene_df_copy


def get_clone_name(
    gene_df: pd.DataFrame,
    data_path: Path=Path("./fine_tune/datasets/genome"),
    allow_accession: bool=False
) -> pd.DataFrame:
    """
    Map accession numbers in gene symbols to clone names using the NCBI Clone DB.

    If a perfect match is found between accession and version, or if the accession is unique,
    the clone name is used. Otherwise, the original symbol is retained. Optionally, unmapped
    accession numbers can be replaced with the Ensembl ID.

    Args:
        gene_df (pd.DataFrame): Input DataFrame with a 'symbol' column containing accession numbers.
        data_path (Path, optional): Directory to store or find the Clone DB table.
        allow_accession (bool, optional): If False, unmapped accession numbers are replaced with Ensembl IDs.

    Returns:
        pd.DataFrame: Copy of the DataFrame with updated 'symbol' and 'CloneName' columns.
    """

    # Description of the table obtained here: https://ftp.ncbi.nlm.nih.gov/repository/clone/reports/000_README_REPORTS.txt

    clone_url = (
        "https://ftp.ncbi.nlm.nih.gov/repository/clone/reports/"
        "Homo_sapiens/clone_acstate_9606.out"
    )
    if not (tmp_clone_path := data_path.joinpath("genome/tmp_clone.txt")).exists():
        download_url(clone_url, tmp_clone_path, print_flag="Clone DB table")
    clone_df = pd.read_csv(tmp_clone_path, sep="\t")

    # Clean up duplicated Accession number.
    #   Stdn: Y: yes standard, A: alias, N: non-standardized/non-aliased name for the clone.
    #   Priority given to standard name: "Y", if no "Y", then "A"
    #   In either case if "Y" or "A" present more than once,
    #   then the accession number is ambiguous and won't be changed.

    print("Clean Clone DB table")
    clone_df = (
        clone_df
        [["Accession", "Stdn", "CloneName"]]
        .groupby("Accession")
        .agg(list)
        .reset_index()
        .loc[lambda df:
            (df["Stdn"].apply(len) == 1) | # if one name, then keep
            (df["Stdn"].apply(lambda x: x.count("Y") == 1)) | # if Y present once, keep
            (df["Stdn"].apply(lambda x: (x.count("A") == 1) & (x.count("Y") == 0))) # if A present once, and not Y, keep
                                                                                    # otherwise ambiguous
        ].explode(["Stdn", "CloneName"])
        .replace({"Stdn": "A"}, "X") # temporarly replace A in X, so that when sorting in decreasing: Y, X, N
        .sort_values("Stdn", ascending=False)
        .drop_duplicates("Accession")
        .replace({"Stdn": "X"}, "A")
        .reset_index(drop=True)
    )

    # Map accession number to clone name:
    #   - if perfect match between accession + version
    #   - if match with accession but accession is unique

    print("Map clone name")
    gene_df_copy = gene_df.copy()
    gene_df_copy = gene_df_copy.set_index("ensembl_id") # move ensembl_id as index
    # Split accession and version
    tmp_split = gene_df_copy["symbol"].str.split(".", expand=True).rename(columns={0: "Accession", 1: "version"})
    tmp_split["version"] = tmp_split["version"].astype(str)
    clone_df[["Accession", "version"]] = clone_df["Accession"].str.split(".", expand=True)
    clone_df["version"] = clone_df["version"].astype(str)

    # Merge on accession and version (when no match, CloneName store NaN)
    merged = (
        tmp_split
        .reset_index(drop=False) # move ensembl_id
        .merge(
            clone_df[["Accession", "version", "CloneName"]],
            on=["Accession", "version"], how="left"
        ).set_index("ensembl_id") # reindex on ensembl_id
    )

    merged["Unique_Accession"] = merged["Accession"].isin(tmp_split[
        tmp_split["Accession"].isin(clone_df["Accession"]) &
        ~tmp_split["Accession"].duplicated(keep=False)
    ]["Accession"])

    # Map CloneName according above logic or retain NaN
    gene_df_copy["CloneName"] = (
        merged.apply(
            lambda row:
                row["CloneName"]
                if pd.notna(row["CloneName"])
                else (
                    clone_df.loc[clone_df["Accession"] == row["Accession"], "CloneName"].iloc[0] # type: ignore
                    if row["Unique_Accession"]
                    else row["CloneName"] # NaN
                ),
            axis=1
        )
    )
    # Unmaped CloneName should fallback to original symbol
    gene_df_copy["symbol"] = gene_df_copy["CloneName"].where(
        gene_df_copy["CloneName"].notna(), # unmaped CloneName should fallback to original symbol
        gene_df_copy["symbol"]
    )

    gene_df_copy = gene_df_copy.reset_index(drop=False) # move back ensembl_id as column
    print(
        "Number of accession number in gene symbol renamed with clone names: "
        f"{(diff := (gene_df_copy['symbol'] != gene_df['symbol']).sum())}")

    # Turn unmaped accession number into ensembl_id
    if not allow_accession:
        gene_df_copy["symbol"] = gene_df_copy["symbol"].where(
            ~gene_df_copy["symbol"].astype(str).str.split(".").apply(lambda x: x[0]).isin(clone_df["Accession"]),
            gene_df_copy["ensembl_id"]
        )
        print("Remaining accession number in gene symbol turned into Ensembl IDs: "
              f"{(gene_df_copy['symbol'] != gene_df['symbol']).sum() - diff}")

    tmp_clone_path.unlink()
    return gene_df_copy


def format_duplicate_symbol(
    gene_df: pd.DataFrame,
    allow_duplicate_symbol: bool=False
) -> pd.DataFrame:
    """
    Handle duplicate or missing gene symbols in a DataFrame.

    For missing gene symbols (NaN), replaces them with the corresponding Ensembl ID.
    For duplicate symbols:
      - If allow_duplicate_symbol is False, replaces all duplicates with their Ensembl IDs.
      - If allow_duplicate_symbol is True, appends a numeric tag to duplicates to make them unique.

    Args:
        gene_df (pd.DataFrame): Input DataFrame with 'symbol' and 'ensembl_id' columns.
        allow_duplicate_symbol (bool, optional): If True, duplicate symbols are indexed with a numeric tag.
            If False, duplicate symbols are replaced with Ensembl IDs.

    Returns:
        pd.DataFrame: DataFrame with updated 'symbol' column, ensuring uniqueness and no missing values.
    """
    gene_df["_symbol"] = gene_df["symbol"].where(
        gene_df["symbol"].notna(),
        gene_df["ensembl_id"]
    )
    nan_count = (gene_df['_symbol'] != gene_df['symbol']).sum()
    if allow_duplicate_symbol:
        def is_tagged(symbol):
            return bool(re.search(r'-\d+$', symbol))

        # Step 1: Sort and prep
        gene_df_copy = gene_df.copy()
        gene_df_copy = (
            gene_df_copy
            .sort_values(["_symbol", "ensembl_id"])
            .reset_index(drop=True)
        )

        ## add a tag to duplicate symbol, and give priority to symbol with already existing tag.
        # Step 2: Group by base symbol (excluding suffixes)
        gene_df_copy["base_symbol"] = gene_df_copy["_symbol"].apply(lambda x: re.sub(r'-\d+$', '', x))
        used_symbols = set(gene_df_copy["_symbol"])  # All existing symbols

        # Step 3: Track renaming per base
        seen = {}
        new_symbols = []

        for _, row in gene_df_copy.iterrows():
            symbol = row["_symbol"]
            base = row["base_symbol"]

            # Leave already-suffixed symbols alone
            if is_tagged(symbol):
                new_symbols.append(symbol)
                continue

            count = seen.get(base, 0)

            if count == 0:
                new_symbol = symbol
            else:
                # Find next available tag not already used
                i = count
                while True:
                    candidate = f"{base}-{i}"
                    if candidate not in used_symbols:
                        break
                    i += 1
                new_symbol = candidate
                used_symbols.add(new_symbol)

            seen[base] = count + 1
            new_symbols.append(new_symbol)

        gene_df_copy["_symbol"] = new_symbols

        # Step 4: Re-integrate
        gene_df = gene_df.set_index("ensembl_id")
        gene_df["_symbol"] = gene_df_copy.set_index("ensembl_id").reindex(gene_df.index)["_symbol"]
        gene_df = gene_df.reset_index()
    else:
        gene_df["_symbol"] = gene_df["_symbol"].where(
            ~gene_df["_symbol"].duplicated(keep=False),
            gene_df["ensembl_id"]
        )
    duplicate_count = (gene_df["symbol"] != gene_df["_symbol"]).sum() - nan_count
    print(
        f"NaNs turned into Ensembl IDs: {nan_count}\n" \
        f"Duplicates {'turned into Ensembl IDs' if not allow_duplicate_symbol else 'indexed with num tag'}: " \
        f"{duplicate_count}"
    )
    gene_df["symbol"] = gene_df["_symbol"]
    return gene_df


def build_ref_genome(
    data_path: Path=Path("./fine_tune/datasets/genome"),
    out_name: Optional[str]=None,
    version: Union[int, str]=32,
    filter_biotype: bool=True,
    use_hgnc: bool=True,
    hgnc_date: Optional[str]="2021-03-01",
    use_clone: bool=True,
    allow_duplicate_symbol: bool=False,
    allow_accession: bool=False
) -> None:
    """
    Download and process a GENCODE GTF annotation file to extract gene IDs and gene names.

    Downloads the GTF file if not present, extracts gene_id and gene_name pairs for all genes
    of selected biotypes, and saves them as a CSV file. Optionally updates gene symbols using
    HGNC and Clone DB mappings, and handles duplicates or accession numbers as specified.

    Args:
        data_path (Path, optional): Directory to store the downloaded and processed files.
        out_name (str, optional): Name for the output CSV file. Defaults to GTF base name.
        version (int or str, optional): GENCODE annotation release version.
        filter_biotype (bool, optional): If True, filter according the following biotype:
            ```
            {"protein_coding", "lncRNA",
            "IG_C_gene", "IG_D_gene", "IG_J_gene", "IG_LV_gene", "IG_V_gene",
            "IG_V_pseudogene", "IG_J_pseudogene", "IG_C_pseudogene",
            "TR_C_gene", "TR_D_gene", "TR_J_gene", "TR_V_gene",
            "TR_V_pseudogene", "TR_J_pseudogene"}
            ```
            And filter out `tag` containing "readthrough_transcript" or "PAR".
        use_hgnc (bool, optional): If True, update gene symbols using HGNC mapping.
        hgnc_date: (str, optional): date of the HGNC achive. Must figure in:
            https://www.genenames.org/download/archive/monthly/tsv/. If None, use_hgnc must be False and conversely.
        use_clone (bool, optional): If True, update gene symbols using Clone DB mapping.
        allow_duplicate_symbol (bool, optional): If True, duplicate symbols are indexed with a numeric tag.
            If False, duplicate symbols are replaced with Ensembl IDs.
        allow_accession (bool, optional): If False, accession numbers are replaced with Ensembl IDs.


    Returns:
        None
    """

    if isinstance(data_path, str):
        data_path = Path(data_path)
    if (use_hgnc and hgnc_date is None):
        raise ValueError(
            "Unexpected `hgnc_date`. `use_hgnc` is set to True, yet `hgnc_date` is None. " \
            "Please set if `use_hgnc` to False or specify a valid date in: " \
            "https://www.genenames.org/download/archive/monthly/tsv/"
        )

    # Set paths
    gtf_url = ("http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
               f"release_{version}/gencode.v{version}.primary_assembly.annotation.gtf.gz")
    tmp_gtf_name = gtf_url.split("/")[-1]
    out_name = out_name or "-".join(tmp_gtf_name.split(".")[:-2])

    # Download GTF if not already downloaded
    out_path = data_path.joinpath(out_name).with_suffix(".csv")
    if out_path.exists():
        print(f"Found local copy of {out_path}")
    else:
        if not (tmp_gtf_path := data_path.joinpath(tmp_gtf_name)).exists():
            download_url(gtf_url, tmp_gtf_path, print_flag=f"Gencode GTF file v{version}")

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

        with gzip.open(tmp_gtf_path, "rt") as gtf:
            print(f"Extracting {tmp_gtf_path}")
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
                    if filter_biotype:
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
                    hgnc_id = attrs.get("hgnc_id", "")
                    if gene_id not in gene_info:
                        gene_info[gene_id] = {"gene_name": gene_name, "hgnc_id": hgnc_id}

        # Retain only allowed genes and create dataframe
        print(f"Retain only allowed genes: {len(gene_allowlist)}")
        filtered_gene_info = {
            gene_id: gene_info.get(gene_id, {"gene_name": "", "hgnc_id": ""}) for gene_id in sorted(gene_allowlist)
        }

        gene_allow_df = (
            pd.DataFrame.from_dict(filtered_gene_info, orient="index")
            .reset_index()
            .rename(columns={"index": "ensembl_id", "gene_name": "symbol"})
        )
        # Update symbol
        if use_hgnc:
            gene_allow_df = get_hgnc_symbol(gene_allow_df, data_path=data_path, hgnc_date=hgnc_date) # type: ignore
        if use_clone:
            gene_allow_df = get_clone_name(gene_allow_df, data_path=data_path, allow_accession=allow_accession)
        # format duplicate: NaN -> ensembl_id, duplicate -> either ensembl_id or indexed with num tag.
        gene_allow_df = format_duplicate_symbol(gene_allow_df, allow_duplicate_symbol=allow_duplicate_symbol)

        gene_allow_df[["ensembl_id", "symbol"]].to_csv(out_path, index=False)
        print(f"Saved {len(gene_allow_df)} gene entries to {out_path}")
        tmp_gtf_path.unlink()

if __name__ == "__main__":
    CLI(build_ref_genome)

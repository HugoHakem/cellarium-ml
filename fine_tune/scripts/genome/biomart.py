import warnings
from functools import reduce
from pathlib import Path
from typing import Optional

import pandas as pd
import scanpy as sc
from jsonargparse import CLI
from pybiomart import Server


def ensembl_to_host(
    ensembl_version: int|str,
    organism: str = "human"
) -> str:
    """
    Get the Ensembl archive host URL for a given Ensembl version and organism.

    Args:
        ensembl_version (int or str): The Ensembl release version.
        organism (str, optional): Organism name (default "human"). Must match an entry in the Ensembl archive table.

    Returns:
        str: The Ensembl archive host URL.

    Raises:
        ValueError: If the organism or Ensembl version is not found in the archive table,
            or if the organism-version pair is not available.

    Example:
        .. code-block:: python

            ensembl_to_host(111, organism="human")
            # 'jan2024.archive.ensembl.org'
    """
    ensembl_version = str(ensembl_version).strip()
    organism = organism.strip().lower()
    url = "https://useast.ensembl.org/info/website/archives/assembly.html"

    # Load all HTML tables from the page
    tables = pd.read_html(url)

    def clean_table(df):
        df = df.rename(columns={"Unnamed: 0": "organism"})
        df = df[df["organism"].notna()]
        df["organism"] = df["organism"].str.strip().str.lower()
        return df.set_index("organism")

    ensembl_archives = reduce(
        lambda left, right: left.join(right, how="outer"),
        map(clean_table, tables)
    )

    if organism not in ensembl_archives.index:
        raise ValueError(
            f"Organism '{organism}' not found in Ensembl archive page:\n{url}"
        )

    ensembl_to_date = (
        ensembl_archives
        .columns
        .str.extract(r"(?P<month>\w+)\s*(?P<year>\d{4})\s*v(?P<version>\d+)", expand=True)
        .dropna()
    )
    ensembl_to_date["archive_id"] = (
        ensembl_to_date["month"].str.lower() + ensembl_to_date["year"].str.lower()
    )
    ensembl_to_date["column"] = ensembl_archives.columns.values
    ensembl_to_date["version"] = ensembl_to_date["version"].astype(str)
    ensembl_to_date.set_index("version", inplace=True)

    if ensembl_version not in ensembl_to_date.index:
        raise ValueError(
            f"Unsupported ensembl_version. You provided ensembl_version: {ensembl_version} "\
            f"but couldn't find it in the table column at {url}"
        )

    col = ensembl_to_date.loc[ensembl_version]["column"]
    if ensembl_archives.loc[organism][col].isna():
        raise ValueError(
            f"The pair (organism={organism}, ensembl_version={ensembl_version}) "
            f"is not available on the archive page: {url}"
        )

    date = ensembl_to_date.loc[ensembl_version]["date"]
    return f"{date}.archive.ensembl.org"


def gencode_to_ensembl(
    gencode_version: int|str,
    organism: str="human"
) -> str:
    """
    Maps a GENCODE version to its corresponding Ensembl release number for a given organism.

    Args:
        gencode_version (int or str): The GENCODE version to look up (e.g., 42).
        organism (str, optional): The organism to query ("human" or "mouse").
            Defaults to "human".

    Returns:
        str: The Ensembl release number corresponding to the given GENCODE version.

    Raises:
        ValueError: If the organism is not supported or the GENCODE version is not found.
        RuntimeError: If there is an error fetching or parsing the GENCODE release table.

    Example:
        .. code-block:: python

            gencode_to_ensembl(42, organism="human")
            # '109'
    """

    gencode_version = str(gencode_version)
    organism = organism.strip().lower()
    if organism not in ["human", "mouse"]:
        raise ValueError(
            f"Invalid organism: {organism}\n"
            "Supported organisms are: 'human', 'mouse'"
        )
    url = f"https://www.gencodegenes.org/{organism}/releases.html"
    try:
        tables = pd.read_html(url)
    except Exception as e:
        raise RuntimeError(f"Failed to fetch or parse HTML from {url}\nError: {e}")

    gencode_release = reduce(
        lambda left, right: left.join(right, how="outer"),
        tables
    )
    gencode_release.columns = [col.strip() for col in gencode_release.columns]
    gencode_release["GENCODE version"] = gencode_release["GENCODE version"].astype(str)
    gencode_release["Ensembl release"] = gencode_release["Ensembl release"].astype(str)

    gencode_release = gencode_release.set_index("GENCODE version")["Ensembl release"]

    if gencode_version not in gencode_release:
        raise ValueError(
            f"GENCODE version {gencode_version} not found in {url}"
        )
    return gencode_release[gencode_version]


def format_organism(organism: str) -> str:
    """
    Maps an organism name to its corresponding Ensembl Biomart dataset prefix.

    Queries the Ensembl Biomart server and matches the input organism against known
    dataset prefixes (e.g., "hsapiens") or display names (e.g., "human").

    Args:
        organism (str): Organism name or dataset prefix (e.g., "human", "mmusculus").

    Returns:
        str: Dataset prefix for use with pybiomart (e.g., "hsapiens").

    Raises:
        ValueError: If the organism name is not recognized.

    Example:
        .. code-block:: python

            format_organism("human")
            # 'hsapiens'

            format_organism("mmusculus")
            # 'mmusculus'
    """
    organism = organism.strip().lower()

    # Connect to Ensembl Biomart
    server = Server(host='http://www.ensembl.org')  # Current release
    mart = server['ENSEMBL_MART_ENSEMBL']
    datasets_df = mart.list_datasets().copy()

    # Extract short name (prefix) and display name
    datasets_df["prefix"] = datasets_df["name"].str.extract(r"(\w+)_gene_ensembl")[0]
    datasets_df["clean_display_name"] = datasets_df["display_name"] \
        .str.extract(r"(.+?) genes", expand=False) \
        .str.strip().str.lower()

    # Build lookup maps
    by_prefix = set(datasets_df["prefix"].dropna())
    by_display = datasets_df.set_index("clean_display_name")["prefix"].dropna()

    if organism in by_prefix:
        return organism
    elif organism in by_display:
        return by_display[organism]
    else:
        raise ValueError(
            f"Invalid organism name: '{organism}'\n"
            "Available names include:\n"
            f"  Prefixes: {sorted(by_prefix)[:5]}...\n"
            f"  Display: {sorted(by_display.index)[:5]}...\n\n"
            "For full list, visit: https://useast.ensembl.org/info/website/archives/assembly.html"
        )


def get_biomart(
    data_path: Path=Path("./fine_tune/datasets/genome"),
    organism: str="human",
    gencode_version: Optional[int]=None,
    ensembl_version: Optional[int|str]=111,
    host: Optional[str]=None,
    attributes: Optional[list[str]]=[
            "ensembl_gene_id",
            "gene_biotype",
            "external_gene_name",
            "external_gene_source",
            "external_synonym"
        ],
    rename_columns: Optional[dict[str, str]]={
            "ensembl_gene_id": "ENSEMBL Gene ID",
            "gene_biotype": "Gene Biotype",
            "external_gene_name": "Gene Symbol",
            "external_gene_source": "Gene Source",
            "external_synonym": "Gene Synonyms"
        },
    print_valid_attributes_on_error: bool = True
) -> pd.DataFrame:
    """
    Downloads and formats Ensembl Biomart gene annotations for a given organism and Ensembl/Gencode version.

    Args:
        data_path (Path): Where to save the gene table (TSV).
        organism (str): Organism name (e.g., "human", "mouse").
        gencode_version (int, optional): GENCODE version to infer Ensembl version.
        ensembl_version (int | str, optional): Ensembl version (if not using GENCODE).
        host (str, optional): Biomart host. If not provided, inferred from Ensembl version.
        attributes (list[str], optional): List of attributes to query.
        rename_columns (bool): Whether to rename columns to more readable names.
        print_valid_attributes_on_error (bool): Whether to print valid attributes if query fails.

    Returns:
        pd.DataFrame: Gene annotation table.
    """
    # Infer Ensembl version if only GENCODE is provided
    if gencode_version is not None and ensembl_version is None:
        ensembl_version = gencode_to_ensembl(gencode_version, organism)

    # Infer host if not given
    if host is None:
        host = ensembl_to_host(ensembl_version, organism)

    # Map to Biomart dataset prefix
    dataset = format_organism(organism)

    try:
        df = sc.queries.biomart_annotations(
            dataset=dataset,
            attrs=attributes,
            host=host,
            use_cache=False
        )
    except Exception as e:
        if print_valid_attributes_on_error:
            print(f"Biomart query failed with error:\n{e}\n")
            print("Listing available attributes for debugging...")
            try:
                server = Server(host=f"http://{host}")
                mart = server['ENSEMBL_MART_ENSEMBL']
                ds = mart.datasets[f"{dataset}_gene_ensembl"]
                attrs_df = ds.attributes
                print(attrs_df[["name", "description"]].to_string(index=False))
            except Exception as inner_e:
                warnings.warn(f"Failed to retrieve valid attributes from host: {inner_e}")
        raise e

    # Optional renaming to user-friendly names
    if rename_columns:
        unmatched_keys = [k for k in rename_columns if k not in df.columns]
        if unmatched_keys:
            warnings.warn(
                "The following columns specified for renaming " \
                f"were not found in the result: {unmatched_keys}"
            )
        df = df.rename(columns={k: v for k, v in rename_columns.items() if k in df.columns})

    # Save to disk
    data_path.mkdir(parents=True, exist_ok=True)
    out_file = data_path / "biomart_gene_table.tsv"
    df.to_csv(out_file, sep="\t", index=False)
    print(f"Biomart annotations saved to: {out_file}")
    return df


if __name__ == "__main__":
    CLI(get_biomart)

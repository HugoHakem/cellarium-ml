import warnings
from functools import reduce
from pathlib import Path
from typing import Any, Optional
from urllib.parse import urlparse

import pandas as pd
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
        organism (str, optional): Organism name (default "human").
            Must match an entry in the Ensembl archive table. Should be a display name.

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
    if not pd.notna(ensembl_archives.loc[organism][col]):
        raise ValueError(
            f"The pair (organism={organism}, ensembl_version={ensembl_version}) "
            f"is not available on the archive page: {url}"
        )

    archive_id = ensembl_to_date.loc[ensembl_version]["archive_id"]
    return f"{archive_id}.archive.ensembl.org"


def gencode_to_ensembl(
    gencode_version: int|str,
    organism: str="human"
) -> str:
    """
    Maps a GENCODE version to its corresponding Ensembl release number for a given organism.

    Args:
        gencode_version (int or str): The GENCODE version to look up (e.g., 42).
        organism (str, optional): The organism to query ("human" or "mouse", in display name).
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
    gencode_release["GENCODE release"] = gencode_release["GENCODE release"].astype(str)
    gencode_release["Ensembl release"] = gencode_release["Ensembl release"].astype(str)

    gencode_release = gencode_release.set_index("GENCODE release")["Ensembl release"]

    if gencode_version not in gencode_release:
        raise ValueError(
            f"GENCODE version {gencode_version} not found in {url}"
        )
    return gencode_release[gencode_version]


def format_organism(organism: str, get_display_name: bool=False, use_cache: bool=False) -> str:
    """
    Maps an organism name to its corresponding Ensembl Biomart dataset prefix or the other way around.

    Args:
        organism (str): Organism name or dataset prefix (e.g., "human", "mmusculus").
        get_display_name (bool): If True, return the display_name, otherwise return the scientific name.
        use_cache (bool): Whether pybiomart should use a cache for requests.
            Will create a .pybiomart.sqlite file in current directory if used.

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
    server = Server(host='http://www.ensembl.org', use_cache=use_cache)  # Latest release
    mart = server['ENSEMBL_MART_ENSEMBL']
    datasets_df = mart.list_datasets().copy()

    # Extract name and display name
    datasets_df["name"] = datasets_df["name"].str.extract(r"(\w+)_gene_ensembl")[0]
    datasets_df["display_name"] = datasets_df["display_name"] \
        .str.extract(r"(.+?) genes", expand=False) \
        .str.strip().str.lower()

    # Build lookup maps
    by_name = datasets_df.set_index("name")["display_name"].dropna()
    by_display = datasets_df.set_index("display_name")["name"].dropna()

    if get_display_name:
        if organism in by_display:
            return organism
        elif organism in by_name:
            return by_name[organism]

    else:
        if organism in by_name:
            return organism
        elif organism in by_display:
            return by_display[organism]

    # If not found, prepare error message
    available_prefixes = sorted(by_name.index)[:5]
    available_display = sorted(by_display.index)[:5]
    raise ValueError(
        f"Invalid organism name: '{organism}'.\n"
        f"Available prefixes (examples): {available_prefixes}...\n"
        f"Available display names (examples): {available_display}...\n"
        "See full list at: https://useast.ensembl.org/info/website/archives/assembly.html"
    )


def _host_parsed(host: str) -> str:
    """
    Parse host url and return dashed separated hostname.

    Args:
        host (str): Biomart host. If not provided, inferred from Ensembl version.

    Returns:
        str: host cleaned up from www. or http:// and from .org

    Example:
        .. code-block:: python

            _host_parsed("http://www.ensembl.org")
            # ensembl
            _host_parsed("http://jan2024.archive.ensembl.org")
            # jan2024-archive-ensembl
    """
    # Ensure URL has scheme for proper parsing
    if not host.startswith(("http://", "https://")):
        host = "https://" + host

    parsed = urlparse(host)
    hostname = parsed.hostname or host  # fallback if parsing fails

    # Remove 'www.' prefix only
    if hostname.startswith("www."):
        hostname = hostname[4:]
    if hostname.endswith(".org"):
        hostname = hostname[:-4]

    # Replace dots with dashes
    host_clean = hostname.replace(".", "-")
    return host_clean


def get_biomart(
    data_path: Path=Path("./fine_tune/datasets/genome"),
    organism: str="human",
    gencode_version: Optional[int]=None,
    ensembl_version: Optional[int|str]=111,
    host: Optional[str]=None,
    attrs: Optional[list[str]|str]=[
            "ensembl_gene_id",
            "gene_biotype",
            "external_gene_name",
            "external_gene_source",
            "external_synonym"
        ],
    filters: dict[str, Any] | None = None,
    rename_columns: Optional[dict[str, str]]={
            "ensembl_gene_id": "ENSEMBL Gene ID",
            "gene_biotype": "Gene Biotype",
            "external_gene_name": "Gene Symbol",
            "external_gene_source": "Gene Source",
            "external_synonym": "Gene Synonyms"
        },
    debug: int=0,
    use_cache: bool=False
) -> pd.DataFrame:
    """
    Downloads and formats Ensembl Biomart gene annotations for a given organism and Ensembl/Gencode version.

    Args:
        data_path (Path): Where to save the gene table (TSV). Can be a filename instead of a directory.
        organism (str): Organism name (e.g., "human", "mouse"). Either display_name or biomart scientific name.
        gencode_version (int, optional): GENCODE version to infer Ensembl version.
        ensembl_version (int | str, optional): Ensembl version (if not using GENCODE).
        host (str, optional): Biomart host. If not provided, inferred from Ensembl version.
        attrs (list[str], str, optional): List of attributes to query.
        rename_columns (dict[str, str], optional): Mapping from original to user-friendly column names.
        debug (int): Level of verbosity when query fails:
            - 0: Silent, just raise raw exception.
            - 1: Show top 10 attributes and filters in the error message.
            - 2: Save full list of attributes/filters to CSV files in `tmp_biomart_output/`.
        use_cache (bool): Whether pybiomart should use a cache for requests.
            Will create a .pybiomart.sqlite file in current directory if used.

    Returns:
        pd.DataFrame: Gene annotation table.
    """
    # Base logic from https://github.com/scverse/scanpy/blob/1.11.2/src/scanpy/queries/_queries.py#L60-L71
    if isinstance(attrs, str):
        attrs = [attrs]
    elif not isinstance(attrs, list):
        raise TypeError(f"attrs must be of type list or str, was {type(attrs)}.")
    try:
        from pybiomart import Server
    except ImportError as e:
        msg = "This method requires the `pybiomart` module to be installed."
        raise ImportError(msg) from e

    # New logic
    ## Infer host using ensembl_version or gencode_version
    if host is None:
        org_display = format_organism(organism=organism, get_display_name=True, use_cache=use_cache)
        org_name = format_organism(organism=organism, get_display_name=False, use_cache=use_cache)
        if gencode_version is not None and ensembl_version is None:
            ensembl_version = gencode_to_ensembl(gencode_version, organism=org_display)
        host = ensembl_to_host(ensembl_version, organism=org_display)

    server = Server(host, use_cache=use_cache)
    dataset = server.marts["ENSEMBL_MART_ENSEMBL"].datasets[f"{org_name}_gene_ensembl"]
    try:
        res = dataset.query(attributes=attrs, filters=filters, use_attr_names=True)
    except Exception as e:
        if debug == 1:
            msg = (
                "Biomart query failed due to invalid `attrs` or `filters`.\n\n"
                "Top available attributes:\n"
                f"{dataset.list_attributes().head(10).to_string(index=False)}\n\n"
                "Top available filters:\n"
                f"{dataset.list_filters().head(10).to_string(index=False)}\n\n"
                "Tip: call `dataset.list_attributes()` and `dataset.list_filters()` to inspect full tables."
            )
        elif debug == 2:
            out_dir = Path("tmp_biomart_output")
            out_dir.mkdir(parents=True, exist_ok=True)
            dataset.list_attributes().to_csv(out_dir / "available_attributes.csv", index=False)
            dataset.list_filters().to_csv(out_dir / "available_filters.csv", index=False)
            msg = (
                "Biomart query failed due to invalid `attrs` or `filters`.\n\n"
                "Full attribute/filter tables saved to:\n"
                f"  - {out_dir / 'available_attributes.csv'}\n"
                f"  - {out_dir / 'available_filters.csv'}"
            )
        else:
            msg = "Biomart query failed. Set `debug=1` or `debug=2` for more info."
        raise Exception(msg) from e

    # ensure every synonym are splitted
    if "external_synonym" in res.columns:
        res["external_synonym"] = (
            res["external_synonym"]
            .apply(lambda x: str(x).split(",") if pd.notnull(x) else [x])
        )
        res = res.explode("external_synonym")
        res["external_synonym"] = res["external_synonym"].str.strip()
        res = res.drop_duplicates()
        res.reset_index(drop=True, inplace=True)


    # Optional renaming to user-friendly names
    if rename_columns:
        unmatched_keys = [k for k in rename_columns if k not in res.columns]
        if unmatched_keys:
            warnings.warn(
                "The following columns specified for renaming " \
                f"were not found in the result: {unmatched_keys}"
            )
        res = res.rename(columns={k: v for k, v in rename_columns.items() if k in res.columns})

    # Save to disk
    if isinstance(data_path, str):
        data_path = Path(data_path)
    data_path.mkdir(parents=True, exist_ok=True)
    if data_path.suffix not in [".csv", ".tsv"]:
        if gencode_version is not None:
            file_name = f"biomart-gencode-v{gencode_version}.csv"
        elif ensembl_version is not None:
            file_name = f"biomart-ensembl-v{ensembl_version}.csv"
        else:
            file_name = f"biomart-host-{_host_parsed(host=host)}.csv"
        data_path = data_path / file_name
    res.to_csv(data_path, index=False)
    print(f"Biomart annotations saved to: {data_path}")
    return res

if __name__ == "__main__":
    CLI(get_biomart)

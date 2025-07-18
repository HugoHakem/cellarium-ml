
import gzip
import logging
import re
import shutil
from collections import defaultdict
from pathlib import Path
from typing import Optional, Union

import numpy as np
import pandas as pd
import requests  # type: ignore
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from intervaltree import Interval, IntervalTree
from tqdm import tqdm

logger = logging.getLogger(__name__)


# === SECTION: Check for tools dependencies ===

class MissingToolError(RuntimeError):
    """
    Exception raised when a required external tool is missing from the system PATH.
    """
    pass

def check_dependencies(tools: str|list[str]):
    """
    Check if required CLI tools are available in the system's PATH.

    Args:
        tools (str or list[str]): Name or list of names of CLI tools to check.

    Raises:
        MissingToolError: If any tool is not found in the system PATH.
    """
    tools = tools if isinstance(tools, list) else [tools]
    missing = [tool for tool in tools if shutil.which(tool) is None]
    if missing:
        raise MissingToolError(
            f"Missing required tools: {', '.join(missing)}.\n"
            f"Please install them and ensure they're in your PATH."
        )


# === SECTION: Download and prepare annotation reference ===


def download_file(url: str, out_path: Path) -> None:
    """
    Download a file from a given URL to a specified local path.

    Args:
        url (str): The URL to download the file from.
        out_path (Path): The local file path where the downloaded file will be saved.

    Returns:
        None
    """
    if out_path.exists():
        logger.info(f"File already exists: {out_path}")
        return
    logger.info(f"Downloading: {url}")
    r = requests.get(url, stream=True)
    with open(out_path, 'wb') as f:
        shutil.copyfileobj(r.raw, f)

def gunzip_file(gz_path: Path, out_path: Optional[Path]=None) -> Path:
    """
    Uncompress a .gz file to a specified output path.

    Args:
        gz_path (Path): Path to the .gz file.
        out_path (Optional[Path]): Path to save the uncompressed file. If None, removes '.gz' extension.

    Returns:
        Path: Path to the uncompressed file.

    Raises:
        ValueError: If out_path is None and gz_path does not end with '.gz'.
    """
    if out_path is None:
        if gz_path.suffix == ".gz":
            out_path = gz_path.with_suffix("")
        else:
            raise ValueError(f"Cannot infer uncompressed file name from {gz_path}")

    if out_path.exists():
        logger.info(f"Uncompressed file already exists: {out_path}")
        return out_path
    logger.info(f"Unzipping {gz_path} → {out_path}")
    with gzip.open(gz_path, 'rb') as f_in, open(out_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    return out_path


def prepare_reference(ref_dir: Path, gencode_version: int) -> tuple[Path, Path]:
    """
    Download and prepare GENCODE reference files (GTF and genome FASTA).

    Args:
        ref_dir (Path): Directory to store reference files and index.
        gencode_version (int): GENCODE release version to use.

    Returns:
        tuple[Path, Path]: Paths to the GTF annotation file and Fasta File.
    """
    ref_dir.mkdir(parents=True, exist_ok=True)
    base_url = f"http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{gencode_version}"

    gtf_gz = ref_dir / "annotation.gtf.gz"
    gtf = gtf_gz.with_suffix("")
    if gtf.exists():
        logger.info(f"Uncompressed file already exists: {gtf}")
    else:
        download_file(f"{base_url}/gencode.v{gencode_version}.primary_assembly.annotation.gtf.gz", gtf_gz)
        gunzip_file(gtf_gz)

    fasta_gz = ref_dir / "genome.fa.gz"
    fasta = fasta_gz.with_suffix("")

    if fasta.exists():
        logger.info(f"Uncompressed file already exists: {fasta}")
    else:
        download_file(f"{base_url}/GRCh38.primary_assembly.genome.fa.gz", fasta_gz)
        gunzip_file(fasta_gz)

    return gtf, fasta


# === SECTION: process gtf file ===


def _find_holes(tree: IntervalTree):
    """
    Find gaps (holes) in an interval tree.

    Args:
        tree (IntervalTree): IntervalTree to search for holes.

    Returns:
        IntervalTree: IntervalTree containing intervals representing holes.
    """
    if not tree:
        return IntervalTree()
    left = tree.begin()
    right = tree.end()
    holes = IntervalTree([Interval(left, right)])

    for iv in tree:
            holes.chop(iv.begin, iv.end)
    return holes


def densify_tree(tree: IntervalTree, split_holes_equally: bool=True) -> IntervalTree:
    """
    Densify an interval tree by filling holes, optionally splitting holes equally.

    Args:
        tree (IntervalTree): IntervalTree to densify.
        split_holes_equally (bool): If True, split holes equally between neighbors.

    Returns:
        IntervalTree: Densified IntervalTree.
    """
    if not tree:
        return IntervalTree()
    holes = sorted(_find_holes(tree))
    densified = tree.copy()
    for hole in holes:
        start, end = hole.begin, hole.end
        # left_offset is here because the interval end is exlusive and in case of odd interval length
        # the middle element must be shared between left and right.
        mid, left_offset = divmod(end + start, 2)
        left_neighbors, right_neighbors = densified[start-1], densified[end]

        for iv in left_neighbors:
            new_iv = Interval(iv.begin, mid + left_offset if split_holes_equally else end, iv.data)
            densified.remove(iv)
            densified.add(new_iv)

        for iv in right_neighbors:
            new_iv = Interval(mid if split_holes_equally else start, iv.end, iv.data)
            densified.remove(iv)
            densified.add(new_iv)
    return densified


def load_gtf_annotations(
    gtf_path: Path,
    densify_gene_interval: bool=True,
    split_holes_equally: bool=True,
) -> dict[str, IntervalTree]:
    """
    Parse a GTF annotation file and extract gene intervals as IntervalTrees.

    Args:
        gtf_path (Path): Path to the GTF annotation file.
        densify_gene_interval (bool): If True, densify gene intervals to guarantee a hit.
        split_holes_equally (bool): If True, split holes equally when densifying.

    Returns:
        dict[str, IntervalTree]: Mapping chromosome names to IntervalTrees of gene intervals.
    """
    logger.info("Parsing GTF file...")
    visited_gene: dict[str, dict[str, dict[str, tuple]]] = defaultdict(
        lambda: defaultdict(
            lambda: defaultdict(tuple)))

    with open(gtf_path) as f:
        for line in tqdm(f, desc="Processing GTF", unit=" lines"):
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            chrom, start, end, info = fields[0], int(fields[3]), int(fields[4]), fields[8]
            attr_pairs = dict(re.findall(r'(\S+) "([^"]+)"', info))
            gene_id = attr_pairs.get("gene_id")
            gene_name = attr_pairs.get("gene_name")
            if gene_id and gene_name:
                if fields[2] in {"gene"}: # to keep the true gene start and end
                    gene_start, gene_end = start, end
                else:
                    gene_start, gene_end = float("inf"), float("-inf") # type: ignore

                if gene_id not in visited_gene[chrom]:
                    visited_gene[chrom][gene_id][gene_name] = (start, end, gene_start, gene_end)
                elif gene_name not in visited_gene[chrom][gene_id]:
                    visited_gene[chrom][gene_id][gene_name] = (start, end)
                else:
                    visited_gene[chrom][gene_id][gene_name] = (
                        min(start, visited_gene[chrom][gene_id][gene_name][0]),
                        max(end, visited_gene[chrom][gene_id][gene_name][1]),
                        min(gene_start, visited_gene[chrom][gene_id][gene_name][2]),
                        max(gene_end, visited_gene[chrom][gene_id][gene_name][3])
                    )
    f.close()
    gene_interval_trees = {chrom: IntervalTree(
        [   # end is exclusive in Interval, so add 1 as gtf start, end annotation are inclusive
            Interval(iv[0], iv[1]+1, (gene_id, gene_name, iv[2], iv[3]))
            for gene_id, gene_name_dict in gene_id_dict.items()
            for gene_name, iv in gene_name_dict.items()
        ])
        for chrom, gene_id_dict in visited_gene.items()}

    if densify_gene_interval:
        return {chr: densify_tree(tree, split_holes_equally) for chr, tree in gene_interval_trees.items()}

    return gene_interval_trees


# === SECTION: Annotation helper functions ===


def _load_targets_csv(targets_path: Path, targets_key_map: Optional[dict] = None) -> pd.DataFrame:
    """
    Load the targets csv file and remap the column as needed.

    Args:
        targets_path (Path): Path to the CSV file containing target sequences.: Union[str, Path],
        targets_key_map (dict, optional): Mapping of required keys ('target_id', 'target_a', 'target_b')
            to actual column names in the CSV if different.
            Example:
              {'target_id': 'id_col', 'target_a': 'a_col', 'target_b': 'b_col'}
    Returns:
        pd.DataFrame: dataframe with at least ['target_id', 'target_a', 'target_b'] as columns.
    """
    df = pd.read_csv(targets_path)

    required_keys = ['target_id', 'target_a', 'target_b']
    missing_keys = [k for k in required_keys if k not in df.columns]

    if missing_keys:
        if targets_key_map is None:
            raise ValueError(
                f"Missing columns in CSV: {missing_keys}. Please provide a `targets_key_map` for these keys."
            )
        if any(k not in targets_key_map for k in missing_keys):
            missing_in_map = [k for k in missing_keys if k not in targets_key_map]
            raise ValueError(
                f"`targets_key_map` is missing mappings for: {missing_in_map}. Provided: {targets_key_map}"
            )
        df = df.rename(columns={v: k for k, v in targets_key_map.items()})
    return df


def targets_to_fasta(
    targets_path: Path,
    output_fasta: Path,
    targets_key_map: Optional[dict] = None
) -> None:
    """
    Convert a CSV file with target sequences into a FASTA file with sequences named {target_id}_A and {target_id}_B.

    Args:
        targets_path (Path): Path to the CSV file containing target sequences.: Union[str, Path],
        output_fasta (str or Path): Path where to save the output FASTA file.
        targets_key_map (dict, optional): Mapping of required keys ('target_id', 'target_a', 'target_b')
            to actual column names in the CSV if different.
            Example:
              {'target_id': 'id_col', 'target_a': 'a_col', 'target_b': 'b_col'}

    Returns:
        None

    Raises:
        ValueError: If required keys are missing from CSV and `targets_key_map` is not provided or incomplete.

    Example:
        .. code-block:: python
        fasta_from_targets("targets.csv", "targets.fa")
        fasta_from_targets("custom.csv", "out.fa", targets_key_map={
            'target_id': 'id_col', 'target_a': 'a_col', 'target_b': 'b_col'
        })
    """
    output_fasta = Path(output_fasta)

    df = _load_targets_csv(targets_path, targets_key_map=targets_key_map)

    records = []
    for _, row in df.iterrows():
        target_id = str(row['target_id'])
        seq_a = Seq(str(row['target_a']))
        seq_b = Seq(str(row['target_b']))
        records.append(SeqRecord(seq=seq_a, id=f"{target_id}_A", description=""))
        records.append(SeqRecord(seq=seq_b, id=f"{target_id}_B", description=""))

    SeqIO.write(records, output_fasta, "fasta")


def find_genes(
    chrom: str, index: Union[int, slice], trees: dict[str, IntervalTree]
) -> list[Interval]:
    """
    Find genes overlapping a given position or interval using an interval tree.

    Args:
        chrom (str): Chromosome name.
        index (int or slice): Genomic position (int) or range (slice).
        trees (dict): Dictionary mapping chromosome names to IntervalTrees.

    Returns:
        list[Interval]: Overlapping intervals, or best guess if no overlap.
            When guessing: if the query falls before the first or after the last interval,
            return the closest one.
    """
    tree = trees.get(chrom)
    if not tree:
        return []

    hits = tree[index]
    if not hits:
        # If the tree is dense, then if there is no hit,
        # either the hit is before the first interval or after the last one.
        first, last = tree.begin(), tree.end()
        try:
            # Case: index is slice or slice-like
            start, stop = index.start, index.stop # type: ignore
            if start is not None:
                if start >= last:
                    return list(tree[last-1])
            if stop is not None:
                if stop <= first:
                    return list(tree[first])
            else:
                return []
        except AttributeError:
            # Case: index is int
            if index < first: # type: ignore
                return list(tree[first]) # type: ignore
            if index >= last: # type: ignore
                return list(tree[last-1]) # type: ignore
            else:
                return []

    return list(hits)


def get_target_length(targets_path: Path, targets_key_map: Optional[dict]=None) -> pd.DataFrame:
    """
    Get the original query lengths for each target sequence in a CSV file.

    Args:
        targets_path (Path): Path to the CSV file containing target sequences.
        targets_key_map (dict, optional): Mapping of required keys ('target_id', 'target_a', 'target_b')
            to actual column names in the CSV if different.
            Example:
              {'target_id': 'id_col', 'target_a': 'a_col', 'target_b': 'b_col'}

    Returns:
        pd.DataFrame: DataFrame with columns 'query_name' and 'original_query_length'.
    """
    sequence_data = _load_targets_csv(targets_path, targets_key_map=targets_key_map)
    sequence_data = sequence_data.melt(
        id_vars="target_id",
        value_vars=["target_a", "target_b"], var_name="target_type", value_name="sequence"
    )
    sequence_data["query_name"] = (
        sequence_data["target_id"].str.cat(
            sequence_data["target_type"].str.extract(r'target_(.*)')[0].str.upper(),
            sep="_"
        )
    )
    sequence_data["original_query_length"] = sequence_data["sequence"].str.len()
    sequence_data = sequence_data[["query_name", "original_query_length"]]
    return sequence_data


def match_annotation_targets(
    annotation_data,
    keep_best_only,
    smallest_interval
):
    """
    Match target_A and target_B pairs by gene and select best or smallest interval matches.

    Args:
        annotation_data (dict): Mapping from query_name to annotation hits.
        keep_best_only (bool): If True, keep only the annotation with the smallest nucleotide mismatch.
        smallest_interval (bool): If True, keep only the annotation with the smallest interval.

    Returns:
        pd.DataFrame: DataFrame of matched annotation results.
    """
    # Match target_A and target_B pairs
    results = []
    keys = sorted(set(k.rsplit("_", 1)[0] for k in annotation_data.keys()))

    for base_key in keys:
        a_hits = annotation_data.get(f"{base_key}_A", [])
        b_hits = annotation_data.get(f"{base_key}_B", [])

        candidates = []
        for g_a in a_hits:
            for g_b in b_hits:
                if g_a["ensembl_id"] == g_b["ensembl_id"]:
                    nm_list: list[float] = [g_a["nm_score"], g_b["nm_score"]] # type: ignore
                    if not np.isnan(nm_list).all():
                        total_nm = np.nansum(nm_list)
                    else:
                        total_nm = np.nan
                    candidates.append({
                        "ensembl_id" : g_a["ensembl_id"],
                        "gene_name": g_a["gene_name"],
                        "total_nm": total_nm,
                        "chr": g_a["chr"],
                        "gene_start": g_a["gene_start"], "gene_end": g_a["gene_end"],
                        "hit_a": g_a["hit"], "hit_b": g_b["hit"],
                    })

        if candidates:
            if keep_best_only:
                best_nm: float = min(c["total_nm"] for c in candidates)
                best_matches = [c for c in candidates if c["total_nm"] == best_nm]
            else:
                best_matches = candidates

            if smallest_interval:
                best_matches = [min(best_matches, key=lambda match: match["gene_end"] - match["gene_start"])]

            for match in best_matches:
                results.append({
                    "target_id": base_key,
                    **match
                })
    return pd.DataFrame(results)


def process_annotated_df(df: pd.DataFrame, targets_path: Path, targets_key_map: Optional[dict]=None) -> pd.DataFrame:
    """
    Post-process the annotated DataFrame to deduplicate records and recover missing target_ids.

    Args:
        df (pd.DataFrame): Annotated DataFrame from annotate_alignments.
        targets_path (Path): Path to the CSV file containing target sequences.
        targets_key_map (dict, optional): Mapping of required keys ('target_id', 'target_a', 'target_b')
            to actual column names in the CSV if different.
            Example:
              {'target_id': 'id_col', 'target_a': 'a_col', 'target_b': 'b_col'}

    Returns:
        pd.DataFrame: DataFrame with all original target_ids, including unmapped ones, and deduplicated records.
    """
    template_df = _load_targets_csv(targets_path, targets_key_map=targets_key_map)

    # Deduplicate full records per target_id
    grouped_rows = []
    for target_id, group in df.groupby("target_id", dropna=False):
        deduped = group.drop_duplicates(subset=[col for col in df.columns if col != "target_id"])
        records = deduped.drop(columns="target_id").to_dict(orient="list")

        # unwrap lists of length 1
        for k, v in records.items():
            if len(v) == 1:
                records[k] = v[0]
        records["target_id"] = target_id
        grouped_rows.append(records)

    final_df = pd.DataFrame(grouped_rows)

    # Outer join to recover all targets (including unmapped)
    merged_df = template_df.merge(final_df, on="target_id", how="left")
    return merged_df

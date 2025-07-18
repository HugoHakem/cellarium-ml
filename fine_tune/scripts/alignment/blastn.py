import logging
import subprocess
import warnings
from collections import defaultdict
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

from .utils import find_genes, get_target_length, load_gtf_annotations, match_annotation_targets, process_annotated_df

logger = logging.getLogger(__name__)

def blastn_build_index(fasta_path: Path, index_path: Path, **kwargs) -> None:
    """
    Build a BLAST nucleotide (blastn) index from a FASTA file if it does not already exist.

    Args:
        fasta_path (Path): Path to the reference FASTA file.
        index_path (Path): Prefix for the blastn index files.
        **kwargs: Additional arguments (currently unused).

    Returns:
        None
    """
    if all(
        index_path.with_suffix(ext).exists() for ext in [
        '.ndb', '.nhr', '.nin', '.njs', '.nog', '.nos', '.not', '.nsq', '.ntf', '.nto'
    ]):
        logger.info(f"blastn index already exists at: {index_path}")
        return
    logger.info("Building blastn index...")
    subprocess.run(
        [
            "makeblastdb",
            "-in", str(fasta_path),
            "-dbtype", "nucl",
            "-parse_seqids",
            "-out", str(index_path)
         ], check=True)


def blastn_align_targets(
    index_path: Path,
    fasta_in: Path,
    output_path: Path,
    n_mismatch: int=1,
    n_thread: int=1,
    **kwargs
) -> None:
    """
    Align target sequences to a reference genome using blastn.

    Args:
        index_path (Path): Prefix for the blastn index files.
        fasta_in (Path): Path to the input FASTA file with target sequences.
        output_path (Path): Path to save the output CSV file.
        n_mismatch (int): here for compatibility with `bowtie_align_targets`. Please discard.
        n_thread (int): Number of threads to use for alignment. Default is 1.
        **kwargs: Additional parameters for blastn (e.g., word_size, num_threads).
            Example includes:
              - word_seize: 12

    Returns:
        None
    """
    logger.info("Aligning targets using blastn...")
    subprocess.run(
        [
            "blastn",
            "-task", "blastn-short",
            "-query", str(fasta_in),
            "-db", str(index_path),
            "-num_threads", str(n_thread),
        ]
        + ([item for option, val in kwargs.items() for item in (f"{option}", f"{val}")] if kwargs else [])
        + [
            "-outfmt", "6", # output format, 6 is tsv
            "-out", str(output_path)
        ],
        check=True
    )
    logger.info(f"csv file written to: {output_path}")


def blastn_process_aligned_seq(
    aligned_targets_path: Path,
    targets_path: Path,
    targets_key_map: Optional[dict]=None,
    n_mismatch: int=1
) -> pd.DataFrame:
    """
    Process BLASTN alignment results and filter alignments based on mismatch criteria.
    Additionally filter out hit where the given chromosomes is not present in both target A and B

    Args:
        aligned_targets_path (Path): Path to the BLASTN output CSV file.
        targets_path (Path): Path to the CSV file containing original target sequences.
        n_mismatch (int): Maximum number of mismatches allowed in alignments. Default is 1.
        targets_key_map (dict, optional): Mapping of required keys ('target_id', 'target_a', 'target_b')
            to actual column names in the CSV if different.
            Example:
              {'target_id': 'id_col', 'target_a': 'a_col', 'target_b': 'b_col'}

    Returns:
        pd.DataFrame: Filtered alignment results with additional columns for mismatch and target type.
    """
    logger.info("Process `blastn` aligned data")
    # Retrieve target original length information
    sequence_data = get_target_length(targets_path=targets_path, targets_key_map=targets_key_map)

    # Filter out aligned_data based n_mismatch
    aligned_data = pd.read_csv(
        aligned_targets_path,
        names=[
            "query_name", "chr", "pident", "length_query",
            "n_mismatch", "gapopen", "query_start", "query_end",
            "alignment_start", "alignment_end", "evalue", "bitscore"
    ], sep="\t")

    aligned_data = aligned_data.merge(sequence_data, how="left", on="query_name")
    aligned_data["total_mismatch"] = (
        aligned_data["n_mismatch"] + aligned_data["original_query_length"] - aligned_data["length_query"]
    )

    aligned_data = aligned_data[aligned_data["total_mismatch"] <= n_mismatch]

    # Filter out hit where chromosomes are not in both targets A and B.
    aligned_data[["gene_id", "target_type"]] = aligned_data["query_name"].str.rsplit("_", n=1, expand=True)

    res = []
    for _, group in aligned_data.groupby(["gene_id"]):
        res.append(
            group[
                (group["chr"]
                .isin(
                    set.intersection(
                        *group.groupby(["target_type"])
                        .agg({"chr": set})["chr"]
                        .to_list()))
        )])
    aligned_data = pd.concat(res, axis=0).reset_index(drop=True)

    return aligned_data


def blastn_annotate_alignments(
    aligned_targets_path: Path,
    targets_path: Path,
    targets_key_map: Optional[dict],
    n_mismatch: int=1,
    gtf_path: Optional[Path]=None,
    gene_interval_trees: Optional[dict]=None,
    densify_gene_interval: bool=True,
    split_holes_equally: bool=True,
    keep_best_only: bool=False,
    smallest_interval: bool=False,
    **kwargs
) -> pd.DataFrame:
    """
    Annotate aligned target sequences with gene information from a GTF file.

    Args:
        aligned_targets_path (Path): Path to the BLASTN output CSV file with alignments.
        targets_path (Path): Path to the CSV file containing original target sequences.
        targets_key_map (dict, optional): Mapping of required keys ('target_id', 'target_a', 'target_b')
            to actual column names in the CSV if different.
            Example:
              {'target_id': 'id_col', 'target_a': 'a_col', 'target_b': 'b_col'}
        n_mismatch (int, optional): Maximum number of mismatches allowed in alignments. Default is 1.
        gtf_path (Path, optional): Path to the GTF annotation file. Used to compute `gene_interval_trees`.
        gene_interval_trees (dict, optional): Precomputed gene interval mapping, as returned by `load_gtf_annotations`.
        densify_gene_interval (bool, optional): If True, densify gene intervals to guarantee a hit. Default is True.
        split_holes_equally (bool, optional): If True, split holes in intervals equally. Default is True.
        keep_best_only (bool, optional): If True, keep only the annotation with the smallest nucleotide mismatch.
        smallest_interval (bool, optional): If True, return only the smallest overlapping interval
            containing both A and B.
        **kwargs: Additional arguments (currently unused).

    Returns:
        pd.DataFrame: Annotated alignment results including Ensembl gene IDs and gene names.

    Raises:
        ValueError: If both `gtf_path` and `gene_interval_trees` are provided, or if neither is provided.
    """
    if (gtf_path is not None) and (gene_interval_trees is not None):
        raise ValueError("Only one of `gtf_path` or `gene_interval_trees` must be provided, not both.")
    if (gtf_path is None) and (gene_interval_trees is None):
        raise ValueError("One of `gtf_path` or `gene_interval_trees` must be provided.")
    if gtf_path is not None:
        gene_interval_trees = load_gtf_annotations(
            gtf_path,
            densify_gene_interval=densify_gene_interval,
            split_holes_equally=split_holes_equally
        )

    aligned_data = blastn_process_aligned_seq(
        aligned_targets_path=aligned_targets_path,
        targets_path=targets_path,
        targets_key_map=targets_key_map,
        n_mismatch=n_mismatch,

    )
    annotation_data = defaultdict(list)
    _warned_densify = True

    logger.info("Retrieve annotation using the GTF file and Aligned file...")
    for _, read in aligned_data.iterrows():
        query_name = read.get("query_name")
        chr = str(read.get("chr"))
        hit_start, hit_end = sorted([
            int(read.get("alignment_start")), # type:ignore
            int(read.get("alignment_end", read.get("alignment_start") + read.get("original_query_length"))) # type:ignore
        ])
        nm_tag = float(read.get("total_mismatch", np.nan))
        gene_ids_names = find_genes(
            chr,
            slice(hit_start, hit_end+1),
            gene_interval_trees # type: ignore
        )
        if not gene_ids_names and _warned_densify:
            warnings.warn(f"No hit has been identified in one `query_name`: {query_name}. You may want to " \
                "densify the gene intervals with `densify_gene_intervals=True` to guarantee a hit.")
            _warned_densify = False
        if gene_ids_names:
            for iv in gene_ids_names:
                gene_id, gene_name, start, end = iv.data
                if gene_id and query_name:
                    annotation_data[query_name].append(
                        {
                            "ensembl_id": gene_id.rsplit(".", 1)[0], # remove version tag if any
                            "gene_name": gene_name,
                            "nm_score": nm_tag,
                            "chr": chr,
                            "hit": hit_start,
                            "gene_start": start,
                            "gene_end": end
                        }
                    )
    annotated_df = match_annotation_targets(annotation_data, keep_best_only, smallest_interval)
    annotated_df = process_annotated_df(annotated_df, targets_path=targets_path, targets_key_map=targets_key_map)
    return annotated_df

import logging
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import pysam

from .utils import find_genes, get_target_length, load_gtf_annotations, match_annotation_targets, process_annotated_df

logger = logging.getLogger(__name__)

def bowtie_build_index(fasta_path: Path, index_path: Path, n_thread: int=1, **kwargs) -> None:
    """
    Build a Bowtie index from a FASTA file if it does not already exist.

    Args:
        fasta_path (Path): Path to the reference FASTA file.
        index_path (Path): Prefix for the Bowtie index files.
        n_thread (int): Number of threads to use for index building. Default is 1.
        **kwargs: Additional arguments for Bowtie (currently unused).

    Returns:
        None
    """
    if all(
        index_path.with_suffix(ext).exists() for ext in [
        '.1.ebwt', '.2.ebwt', '.3.ebwt', '.4.ebwt', '.rev.1.ebwt', '.rev.2.ebwt'
    ]):
        logger.info(f"Bowtie index already exists at: {index_path}")
        return
    logger.info("Building bowtie index...")
    index_path.parent.mkdir(parents=True, exist_ok=True)
    subprocess.run(["bowtie-build", "--threads", str(n_thread), str(fasta_path), str(index_path)], check=True)


def bowtie_align_targets(
    index_path: Path,
    fasta_in: Path,
    output_path: Path,
    n_mismatch: int=1,
    n_thread: int=1,
    **kwargs
) -> None:
    """
    Align target sequences to a reference genome using Bowtie and output SAM format.

    Args:
        index_path (Path): Prefix for the Bowtie index files.
        fasta_in (Path): Path to the input FASTA file with target sequences.
        output_path (Path): Path to save the output SAM file.
        n_mismatch (int): Maximum number of mismatches allowed in alignments. Default is 1.
        n_thread (int): Number of threads to use for alignment. Default is 1.
        **kwargs: Additional parameters for Bowtie (e.g., custom flags).

    Returns:
        None

    Raises:
        RuntimeError: If Bowtie or Samtools commands fail.
    """
    logger.info("Aligning targets using bowtie + samtools...")

    # Bowtie command (with --sam) piped into samtools view -h
    bowtie_cmd = (
        ["bowtie"]
        + ([item for option, val in kwargs.items() for item in (f"{option}", f"{val}")] if kwargs else [])
        + [
            "-v", str(n_mismatch),
            "--threads", str(n_thread),
            "-a",
            "--best",
            "--strata",
            "-x", str(index_path),
            "-f", str(fasta_in),
            "--sam"
        ]
    )

    samtools_cmd = ["samtools", "view", "-h", "-"]

    with open(output_path, "w") as f_out:
        p1 = subprocess.Popen(bowtie_cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(samtools_cmd, stdin=p1.stdout, stdout=f_out)
        # allow p1 to receive SIGPIPE if p2 exits
        p1.stdout.close() # type: ignore
        p2.communicate()

        if p1.wait() != 0:
            raise RuntimeError("Bowtie command failed")
        if p2.returncode != 0:
            raise RuntimeError("Samtools conversion failed")

    logger.info(f"SAM file written to: {output_path}")


def bowtie_annotate_alignments(
    aligned_targets_path: Path,
    targets_path: Path,
    targets_key_map: Optional[dict],
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
        aligned_targets_path (Path): Path to the SAM file with alignments.
        targets_path (Path): Path to the CSV file containing original target sequences.
        targets_key_map (dict, optional): Mapping of required keys ('target_id', 'target_a', 'target_b')
            to actual column names in the CSV if different.
            Example:
              {'target_id': 'id_col', 'target_a': 'a_col', 'target_b': 'b_col'}
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
    if gene_interval_trees is None:
        gene_interval_trees = load_gtf_annotations(
            gtf_path, # type: ignore
            densify_gene_interval=densify_gene_interval,
            split_holes_equally=split_holes_equally
        )

    samfile = pysam.AlignmentFile(str(aligned_targets_path), "r")
    sequence_data = get_target_length(targets_path=targets_path)
    sequence_data.set_index("query_name", inplace=True)
    annotation_data = defaultdict(list)
    _warned_densify = True

    logger.info("Retrieve annotation using the GTF file and SAM file...")
    for read in samfile:
        if read.is_unmapped:
            continue
        query_name = read.query_name
        chr = str(read.reference_name)
        hit = read.reference_start
        nm_tag = float(read.get_tag("NM")) if read.has_tag("NM") else np.nan
        # bowtie return hit at the forward strand, so we check the `sequence_length` nucleotide ahead as well
        # (the whole match)
        gene_ids_names = find_genes(
            chr,
            slice(hit, hit+sequence_data.loc[query_name]["original_query_length"]),
            gene_interval_trees
        )
        if not gene_ids_names and _warned_densify and (chr in gene_interval_trees):
            logger.warning(f"No hit has been identified in one `query_name`: {query_name}. You may want to " \
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
                            "hit": hit,
                            "gene_start": start,
                            "gene_end": end
                        }
                    )
    annotated_df = match_annotation_targets(annotation_data, keep_best_only, smallest_interval)
    annotated_df = process_annotated_df(annotated_df, targets_path=targets_path, targets_key_map=targets_key_map)
    return annotated_df

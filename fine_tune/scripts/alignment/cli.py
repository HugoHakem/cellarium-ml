
import logging
import shutil
from pathlib import Path
from typing import Callable, Optional, Union

from jsonargparse import CLI

from fine_tune.scripts.alignment.blastn import blastn_align_targets, blastn_annotate_alignments, blastn_build_index
from fine_tune.scripts.alignment.bowtie import bowtie_align_targets, bowtie_annotate_alignments, bowtie_build_index
from fine_tune.scripts.alignment.utils import (
    check_dependencies,
    load_gtf_annotations,
    prepare_reference,
    targets_to_fasta,
)

logger = logging.getLogger(__name__)

def _setup_logger():
    logging.basicConfig(
        level=logging.INFO,
        format="[%(asctime)s] %(levelname)s %(name)s - %(message)s",
        datefmt="%H:%M:%S"
    )


def main(
    targets_paths: Union[Path,list[Path]],
    targets_key_maps: Union[Optional[dict], list[Optional[dict]]]=None,
    out_dir: Path=Path("output"),
    gencode_version: int=45,
    tool: str="bowtie",
    aligner_kwargs: Optional[dict]=None,
    n_mismatch: int=1,
    n_thread: int=1,
    densify_gene_interval: bool=True,
    split_holes_equally: bool=True,
    keep_best_only: bool=False,
    smallest_interval: bool=False,
    cleanup: bool=False
):
    """
    Run the full workflow for aligning target sequences to a reference genome, annotating, and saving results.

    This function processes one or more CSV files containing two target sequences per gene,
    builds a genome index, converts targets to FASTA, aligns them using Bowtie or BLASTN,
    annotates the alignments with gene information from a GTF file, deduplicates and merges results,
    and writes the final annotated CSVs to disk.

    Args:
        targets_paths (Path or list[Path]): Path(s) to input CSV file(s) with target sequences.
        targets_key_maps (dict or list of dict, optional): Mapping(s) of required keys
            ('target_id', 'target_a', 'target_b') to actual column names in the CSV(s) if different.
            Example: {'target_id': 'id_col', 'target_a': 'a_col', 'target_b': 'b_col'}
        out_dir (Path, optional): Output directory for results and intermediate files. Defaults to "output".
        gencode_version (int, optional): GENCODE release version to use for reference files. Defaults to 45.
        tool (str, optional): Alignment tool to use ("bowtie" or "blastn"). Defaults to "bowtie".
        aligner_kwargs (dict, optional): Additional keyword arguments for the aligner.
        n_mismatch (int, optional): Maximum number of mismatches allowed in alignments. Defaults to 1.
        n_thread (int, optional): Number of threads to use for index building and alignment. Defaults to 1.
        densify_gene_interval (bool, optional): If True, densify gene intervals to guarantee a hit. Defaults to True.
        split_holes_equally (bool, optional): If True, split holes in intervals equally. Defaults to True.
        keep_best_only (bool, optional): If True, keep only the annotation with the smallest nucleotide mismatch.
        smallest_interval (bool, optional): If True, return only the smallest overlapping interval
            containing both A and B.
        cleanup (bool, optional): If True, remove intermediate files after processing. Defaults to False.

    Returns:
        None
    """
    out_dir = Path(out_dir)
    tmp_dir = out_dir / "tmp"
    ref_dir = tmp_dir / f"reference/gencode_v{gencode_version}"
    if tool == "bowtie":
        check_dependencies(["bowtie", "samtools"])
        index_builder: Callable = bowtie_build_index
        aligner: Callable = bowtie_align_targets
        annotator: Callable = bowtie_annotate_alignments
        index_path = ref_dir / "bowtie/genome_index"
        _aligned_targets_path_ext = ".sam"
    elif tool == "blastn":
        check_dependencies(["blastn", "makeblastdb"])
        index_builder = blastn_build_index
        aligner = blastn_align_targets
        annotator = blastn_annotate_alignments
        index_path = ref_dir / "blastn/genome_index"
        _aligned_targets_path_ext = ".tsv"
    else:
        raise ValueError(f"Unsupported tool: {tool}")

    gtf, fasta = prepare_reference(ref_dir, gencode_version)

    index_builder(
        fasta,
        index_path,
        n_thread=n_thread
    )

    gene_interval_trees = load_gtf_annotations(
        gtf,
        densify_gene_interval=densify_gene_interval,
        split_holes_equally=split_holes_equally
    )

    targets_paths = targets_paths if isinstance(targets_paths, list) else [targets_paths]
    targets_key_maps = targets_key_maps if isinstance(targets_key_maps, list) else [targets_key_maps]

    for targets_path, targets_key_map in zip(targets_paths, targets_key_maps):
        targets_path = Path(targets_path)
        file_name = targets_path.stem

        fasta_targets_path = tmp_dir / f"{file_name}_targets.fa"
        targets_to_fasta(targets_path, fasta_targets_path, targets_key_map=targets_key_map)

        aligned_targets_path = tmp_dir.joinpath(f"{file_name}_{tool}_aligned").with_suffix(_aligned_targets_path_ext)
        if aligner_kwargs:
            aligner(
                index_path=index_path,
                fasta_in=fasta_targets_path,
                output_path=aligned_targets_path,
                n_mismatch=n_mismatch,
                n_thread=n_thread,
                **aligner_kwargs
            )
        else:
            aligner(
                index_path=index_path,
                fasta_in=fasta_targets_path,
                output_path=aligned_targets_path,
                n_mismatch=n_mismatch,
                n_thread=n_thread
            )

        annotated_df = annotator(
            aligned_targets_path=aligned_targets_path,
            targets_path=targets_path,
            targets_key_map=targets_key_map,
            n_mismatch=n_mismatch,
            gene_interval_trees=gene_interval_trees,
            densify_gene_interval=densify_gene_interval,
            split_holes_equally=split_holes_equally,
            keep_best_only=keep_best_only,
            smallest_interval=smallest_interval
        )

        output_csv = out_dir / f"{file_name}_aligned_{tool}_gencode_v{gencode_version}.csv"
        annotated_df.to_csv(output_csv, index=False)
        logger.info(f"Saved annotated results to {output_csv}")

    if cleanup:
        logger.info("Cleaning up intermediate files...")
        try:
            shutil.rmtree(tmp_dir)
            logger.info(f"Deleted {tmp_dir}")
        except Exception as e:
            logger.error(f"Could not delete {tmp_dir}: {e}")


if __name__ == "__main__":
    _setup_logger()
    CLI(main)

import gzip
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Optional, Union

import pandas as pd
import pysam
import requests  # type: ignore
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from jsonargparse import CLI

# === SECTION: Check for tools dependencies ===

class MissingToolError(RuntimeError):
    """
    Exception raised when a required external tool is missing from the system PATH.
    """
    pass

def _check_dependencies(tools: str|list[str]):
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
        print(f"File already exists: {out_path}")
        return
    print(f"Downloading: {url}")
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
        print(f"Uncompressed file already exists: {out_path}")
        return out_path
    print(f"Unzipping {gz_path} → {out_path}")
    with gzip.open(gz_path, 'rb') as f_in, open(out_path, 'wb') as f_out:
        shutil.copyfileobj(f_in, f_out)
    return out_path


def build_bowtie_index(fasta_path: Path, index_prefix: Path) -> None:
    """
    Build a Bowtie index from a FASTA file if it does not already exist.

    Args:
        fasta_path (Path): Path to the reference FASTA file.
        index_prefix (Path): Prefix for the Bowtie index files.

    Returns:
        None
    """
    if all(
        Path(f"{index_prefix}.{ext}").exists() for ext in [
        '1.ebwt', '2.ebwt', '3.ebwt', '4.ebwt', 'rev.1.ebwt', 'rev.2.ebwt'
    ]):
        print(f"Bowtie index already exists at: {index_prefix}")
        return
    print("Building bowtie index...")
    subprocess.run(["bowtie-build", str(fasta_path), str(index_prefix)], check=True)


def prepare_reference(ref_dir: Path, gencode_version: int) -> tuple[Path, Path]:
    """
    Download and prepare GENCODE reference files (GTF and genome FASTA), and build Bowtie index.

    Args:
        ref_dir (Path): Directory to store reference files and index.
        gencode_version (int): GENCODE release version to use.

    Returns:
        tuple[Path, Path]: Paths to the GTF annotation file and Bowtie index prefix.
    """
    ref_dir.mkdir(parents=True, exist_ok=True)
    base_url = f"http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{gencode_version}"

    gtf_gz = ref_dir / "annotation.gtf.gz"
    gtf = gtf_gz.with_suffix("")
    if gtf.exists():
        print(f"Uncompressed file already exists: {gtf}")
    else:
        download_file(f"{base_url}/gencode.v{gencode_version}.primary_assembly.annotation.gtf.gz", gtf_gz)
        gunzip_file(gtf_gz)

    fasta_gz = ref_dir / "genome.fa.gz"
    fasta = fasta_gz.with_suffix("")

    if fasta.exists():
        print(f"Uncompressed file already exists: {fasta}")
    else:
        download_file(f"{base_url}/gencode.v{gencode_version}.primary_assembly.genome.fa.gz", fasta_gz)
        gunzip_file(fasta_gz)

    index_prefix = ref_dir / "genome_index"
    build_bowtie_index(fasta, index_prefix)

    return gtf, index_prefix


# === SECTION: nnotate csv file ===


def csv_to_fasta(
    csv_path: Union[str, Path],
    output_fasta: Union[str, Path],
    key_map: Optional[dict] = None
) -> None:
    """
    Convert a CSV file with target sequences into a FASTA file with sequences
    named {gene_id}_A and {gene_id}_B.

    Args:
        csv_path (str or Path): Path to the input CSV file.
        output_fasta (str or Path): Path where to save the output FASTA file.
        key_map (dict, optional): Mapping of required keys ('gene_id', 'target_a', 'target_b')
            to actual column names in the CSV if different. Example:
            {'gene_id': 'id_col', 'target_a': 'a_col', 'target_b': 'b_col'}

    Returns:
        None

    Raises:
        ValueError: If required keys are missing from CSV and `key_map` is not provided or incomplete.

    Example:
        .. code-block:: python
        fasta_from_targets("targets.csv", "targets.fa")
        fasta_from_targets("custom.csv", "out.fa", key_map={
            'gene_id': 'id_col', 'target_a': 'a_col', 'target_b': 'b_col'
        })
    """
    csv_path = Path(csv_path)
    output_fasta = Path(output_fasta)

    df = pd.read_csv(csv_path)

    required_keys = ['gene_id', 'target_a', 'target_b']
    missing_keys = [k for k in required_keys if k not in df.columns]

    if missing_keys:
        if key_map is None:
            raise ValueError(
                f"Missing columns in CSV: {missing_keys}. Please provide a `key_map` for these keys."
            )
        if any(k not in key_map for k in missing_keys):
            missing_in_map = [k for k in missing_keys if k not in key_map]
            raise ValueError(
                f"`key_map` is missing mappings for: {missing_in_map}. Provided: {key_map}"
            )
        df = df.rename(columns={v: k for k, v in key_map.items()})

    records = []
    for _, row in df.iterrows():
        gene_id = str(row['gene_id'])
        seq_a = Seq(str(row['target_a']))
        seq_b = Seq(str(row['target_b']))
        records.append(SeqRecord(seq=seq_a, id=f"{gene_id}_A", description=""))
        records.append(SeqRecord(seq=seq_b, id=f"{gene_id}_B", description=""))

    SeqIO.write(records, output_fasta, "fasta")


def align_targets(index_prefix: Path, fasta_in: Path, sam_out: Path, n_mismatch: int=1) -> None:
    """
    Align target sequences to a reference genome using Bowtie and output SAM format.

    Args:
        index_prefix (Path): Prefix for the Bowtie index files.
        fasta_in (Path): Path to the input FASTA file with target sequences.
        sam_out (Path): Path to save the output SAM file.
        n_mismatch (int): Maximum number of mismatches allowed in alignments.

    Returns:
        None

    Raises:
        RuntimeError: If Bowtie or Samtools commands fail.
    """
    print("Aligning targets using bowtie + samtools...")

    # Bowtie command (with --sam) piped into samtools view -h
    bowtie_cmd = [
        "bowtie", "-v", str(n_mismatch), "-a", "--best", "--strata",
        "-x", str(index_prefix),
        "-f", str(fasta_in),
        "--sam"
    ]
    samtools_cmd = ["samtools", "view", "-h", "-"]

    with open(sam_out, "w") as f_out:
        p1 = subprocess.Popen(bowtie_cmd, stdout=subprocess.PIPE)
        p2 = subprocess.Popen(samtools_cmd, stdin=p1.stdout, stdout=f_out)
        # allow p1 to receive SIGPIPE if p2 exits
        p1.stdout.close() # type: ignore
        p2.communicate()

        if p1.wait() != 0:
            raise RuntimeError("Bowtie command failed")
        if p2.returncode != 0:
            raise RuntimeError("Samtools conversion failed")

    print(f"SAM file written to: {sam_out}")


def load_gtf_annotations(gtf_path: Path) -> dict:
    """
    Parse a GTF annotation file and extract gene intervals.

    Args:
        gtf_path (Path): Path to the GTF annotation file.

    Returns:
        dict: Dictionary mapping chromosome names to lists of gene intervals (start, end, gene_id, gene_name).
    """
    print("Parsing GTF file...")
    gene_intervals = defaultdict(list)

    with open(gtf_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split('\t')
            if fields[2] not in {"gene", "transcript"}:
                continue
            chrom, start, end, info = fields[0], int(fields[3]), int(fields[4]), fields[8]
            attr_pairs = dict(re.findall(r'(\S+) "([^"]+)"', info))
            gene_id = attr_pairs.get("gene_id")
            gene_name = attr_pairs.get("gene_name")
            if gene_id and gene_name:
                gene_intervals[chrom].append((start, end, gene_id, gene_name))

    return gene_intervals


def find_gene(chrom, pos, gene_intervals) -> tuple[str|None, str|None]:
    """
    Find the gene overlapping a given chromosome and position.

    Args:
        chrom (str): Chromosome name.
        pos (int): Genomic position (1-based).
        gene_intervals (dict): Dictionary of gene intervals as returned by load_gtf_annotations.

    Returns:
        tuple[str|None, str|None]: (gene_id, gene_name) if found, otherwise (None, None).
    """
    for start, end, gene_id, gene_name in gene_intervals.get(chrom, []):
        if start <= pos <= end:
            return gene_id, gene_name
    return None, None


def annotate_alignments(sam_path: Path, gtf_path: Path) -> pd.DataFrame:
    """
    Annotate aligned target sequences with gene information from a GTF file.

    Args:
        sam_path (Path): Path to the SAM file with alignments.
        gtf_path (Path): Path to the GTF annotation file.

    Returns:
        pd.DataFrame: DataFrame with annotated alignment results, including gene IDs and names.
    """
    gene_intervals = load_gtf_annotations(gtf_path)
    samfile = pysam.AlignmentFile(str(sam_path), "r")
    alignment_data = defaultdict(list)

    for read in samfile:
        if read.is_unmapped:
            continue
        query_name = read.query_name
        chrom = read.reference_name
        pos = read.reference_start
        nm_tag = read.get_tag("NM") if read.has_tag("NM") else 0
        gene_id, gene_name = find_gene(chrom, pos, gene_intervals)
        if gene_id and query_name:
            alignment_data[query_name].append(
                {
                    "ensembl_id": gene_id.rsplit(".", 1)[0], # remove version tag if any
                    "gene_name": gene_name,
                    "nm_score": nm_tag,
                    "chrom": chrom,
                    "pos": pos
                }
            )
    # Match target_A and target_B pairs
    results = []
    keys = sorted(set(k.rsplit("_", 1)[0] for k in alignment_data.keys()))

    for base_key in keys:
        a_hits = alignment_data.get(f"{base_key}_A", [])
        b_hits = alignment_data.get(f"{base_key}_B", [])

        candidates = []
        for g_a in a_hits:
            for g_b in b_hits:
                if g_a["ensembl_id"] == g_b["ensembl_id"]:
                    total_nm = g_a["nm_score"] + g_b["nm_score"] # type: ignore
                    candidates.append({
                        "ensembl_id" : g_a["ensembl_id"],
                        "gene_name": g_a["gene_name"],
                        "total_nm": total_nm,
                        "chrom_a": g_a["chrom"], "pos_a": g_a["pos"],
                        "chrom_b": g_b["chrom"], "pos_b": g_b["pos"],
                    })
        if candidates:
            best_nm = min(c["total_nm"] for c in candidates) # type: ignore
            best_matches = [c for c in candidates if c["total_nm"] == best_nm]
            for match in best_matches:
                results.append({
                    "gene_id": base_key,
                    **match
                })

    return pd.DataFrame(results)


def process_annotated_df(df: pd.DataFrame, template_csv_path: Path) -> pd.DataFrame:
    """
    Post-process the annotated DataFrame to deduplicate records and recover missing gene_ids.

    Args:
        df (pd.DataFrame): Annotated DataFrame from annotate_alignments.
        template_csv_path (Path): Path to the original CSV file with all target gene_ids.

    Returns:
        pd.DataFrame: DataFrame with all original gene_ids, including unmapped ones, and deduplicated records.
    """
    template_df = pd.read_csv(template_csv_path)

    # Deduplicate full records per target_id
    grouped_rows = []
    for gene_id, group in df.groupby("gene_id", dropna=False):
        deduped = group.drop_duplicates(subset=[col for col in df.columns if col != "gene_id"])
        records = deduped.drop(columns="gene_id").to_dict(orient="list")

        # unwrap lists of length 1
        for k, v in records.items():
            if len(v) == 1:
                records[k] = v[0]
        records["gene_id"] = gene_id
        grouped_rows.append(records)

    final_df = pd.DataFrame(grouped_rows)

    # Outer join to recover all targets (including unmapped)
    merged_df = template_df.merge(final_df, on="gene_id", how="left")
    return merged_df


# === SECTION: main logic ===


def main(
    csv_paths: Path|list[Path],
    out_root: Path=Path("output"),
    gencode_version: int=45,
    n_mismatch: int=1,
    cleanup: bool=False,
):
    """
    Run the full workflow for aligning target sequences to a reference genome, annotating, and saving results.

    This function processes one or more CSV files containing target sequences, aligns them to a reference genome
    using Bowtie, annotates the alignments with gene information from a GTF file, deduplicates and merges results,
    and writes the final annotated CSVs to disk. Optionally, it cleans up intermediate files.

    Args:
        csv_paths (Path or list[Path]): Path(s) to input CSV file(s) with target sequences.
        out_root (Path, optional): Output directory for results and intermediate files. Defaults to "output".
        gencode_version (int, optional): GENCODE release version to use for reference files. Defaults to 45.
        n_mismatch (int, optional): Maximum number of mismatches allowed in alignments. Defaults to 1.
        cleanup (bool, optional): If True, remove intermediate files after processing. Defaults to False.

    Returns:
        None
    """
    out_root = Path(out_root)
    tmp_dir = out_root / "tmp"
    ref_dir = tmp_dir / "reference"
    gtf, index_prefix = prepare_reference(ref_dir, gencode_version)

    csv_paths = csv_paths if isinstance(csv_paths, list) else [csv_paths]

    for csv_path in csv_paths:
        csv_path = Path(csv_path)
        sample_name = csv_path.stem

        fasta_targets = tmp_dir / f"{sample_name}_targets.fa"
        csv_to_fasta(csv_path, fasta_targets)

        sam_out = tmp_dir / f"{sample_name}_aligned.sam"
        align_targets(index_prefix, fasta_targets, sam_out, n_mismatch=n_mismatch)

        annotated_df = annotate_alignments(sam_out, gtf)
        annotated_df = process_annotated_df(annotated_df, template_csv_path=csv_path)

        output_csv = out_root / f"{sample_name}_aligned_gencode_v{gencode_version}.csv"
        annotated_df.to_csv(output_csv, index=False)
        print(f"Saved annotated results to {output_csv}")

    if cleanup:
        print("Cleaning up intermediate files...")
        try:
            shutil.rmtree(tmp_dir)
            print(f"Deleted {tmp_dir}")
        except Exception as e:
            print(f"Could not delete {tmp_dir}: {e}")


if __name__ == "__main__":
    _check_dependencies(["bowtie", "samtools"])
    CLI(main)

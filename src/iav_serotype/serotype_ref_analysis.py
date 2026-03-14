import os
import argparse
import subprocess
import logging
from pathlib import Path
import polars as pl
import tempfile
import pysam

logger = logging.getLogger("iavs_logger")


def filter_bam_by_quality(input_bam: str, output_bam: str) -> None:
    """
    Filter BAM file using the same quality criteria as iav_serotype:
    - Alignment length >= 100
    - Alignment proportion (aligned length / read length) >= 0.9
    - Alignment accuracy (1 - edit distance / aligned length) >= 0.8
    """
    with pysam.AlignmentFile(input_bam, 'rb') as in_bam, \
         pysam.AlignmentFile(output_bam, 'wb', template=in_bam) as out_bam:

        for record in in_bam:
            read_length = record.infer_read_length()
            align_length = record.query_alignment_length

            # Skip if read length or alignment length is invalid
            if read_length == 0 or align_length == 0:
                continue

            # Calculate alignment proportion
            try:
                align_prop = align_length / read_length
            except ZeroDivisionError:
                align_prop = 0.0

            # Calculate alignment accuracy using NM tag (edit distance)
            try:
                nm_tag = dict(record.get_tags(with_value_type=False)).get('NM', 0)
                align_acc = (align_length - nm_tag) / align_length
            except ZeroDivisionError:
                align_acc = 0.0

            # Apply quality filters matching iav_serotype's criteria
            if align_length >= 100 and align_prop >= 0.9 and align_acc >= 0.8:
                out_bam.write(record)

    logger.info(f"Filtered BAM saved to: {output_bam}")


def run_samtools_coverage(bam_path: str, output_tsv: str) -> None:
    """Execute samtools coverage to generate alignment coverage statistics"""
    cmd = [
        "samtools", "coverage",
        "-o", output_tsv,
        bam_path
    ]
    try:
        subprocess.run(
            cmd,
            check=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True
        )
        logger.info(f"Generated coverage file: {output_tsv}")
    except subprocess.CalledProcessError as e:
        logger.error(f"samtools coverage failed: {e.stderr}")
        raise


def filter_coverage(coverage_tsv: str, filtered_tsv: str) -> None:
    """Filter coverage data to retain entries with coverage ≥ 1"""
    coverage_headers = [
        "rname", "startpos", "endpos", "numreads",
        "covbases", "coverage", "meandepth", "meanbaseq", "meanmapq"
    ]

    with open(coverage_tsv, 'r') as f:
        lines = [line.strip() for line in f if not line.startswith('#') and line.strip()]

    with tempfile.NamedTemporaryFile(mode='w', delete=False) as tmp:
        tmp.write('\t'.join(coverage_headers) + '\n')
        tmp.write('\n'.join(lines) + '\n')
        tmp_path = tmp.name

    df = pl.read_csv(
        tmp_path,
        separator='\t',
        has_header=True,
        schema_overrides={"meanmapq": pl.Float64}  # Explicitly set meanmapq to float
    )
    os.unlink(tmp_path)

    filtered_df = df.filter(pl.col("coverage") >= 1.0)
    filtered_df.write_csv(filtered_tsv, separator='\t', include_header=True)
    logger.info(f"Filtered coverage file saved to: {filtered_tsv}")


def add_reference_info(filtered_coverage: str, db_info: str, output_tsv: str) -> None:
    """Merge coverage data with reference genome metadata"""
    coverage_df = pl.read_csv(
        filtered_coverage,
        separator='\t',
        has_header=True
    )

    db_df = pl.read_csv(
        db_info,
        separator='\t',
        has_header=True,
        schema_overrides={
            "accession": pl.Utf8,
            "serotype": pl.Utf8,
            "segment": pl.Utf8,
            "Organism_Name": pl.Utf8,
            "Host": pl.Utf8,
            "Collection_Date": pl.Utf8
        }
    )

    merged_df = coverage_df.join(
        db_df,
        left_on="rname",
        right_on="accession",
        how="inner"
    )

    if "accession" in merged_df.columns:
        merged_df = merged_df.drop("accession")

    if merged_df.height == 0:
        logger.warning("No matching reference genome information found in database")

    merged_df.write_csv(output_tsv, separator='\t', include_header=True)
    logger.info(f"Annotated coverage file saved to: {output_tsv}")


def generate_consensus(bam_path: str, ref_fasta: str, filtered_coverage: str,
                      output_vcf: str, output_fasta: str, cov_thresh: int = 10) -> None:
    """Generate consensus sequence with low-coverage regions masked as 'N'"""
    coverage_df = pl.read_csv(
        filtered_coverage,
        separator='\t',
        has_header=True
    )

    if coverage_df.height == 0:
        logger.warning("No references with sufficient coverage. Skipping consensus generation.")
        return

    ref_names = coverage_df["rname"].to_list()
    logger.info(f"Generating consensus for {len(ref_names)} references with coverage ≥ 1")

    # Create temporary reference file with high-coverage sequences
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp_ref:
        tmp_ref_path = tmp_ref.name
        cmd = ["samtools", "faidx", ref_fasta] + ref_names
        try:
            result = subprocess.run(
                cmd,
                check=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )
            tmp_ref.write(result.stdout)
        except subprocess.CalledProcessError as e:
            logger.error(f"Failed to extract reference sequences: {e.stderr}")
            os.unlink(tmp_ref_path)
            raise

    # Create temporary directory for intermediate files
    with tempfile.TemporaryDirectory() as tmp_dir:
        # Generate variants and initial consensus
        vcf_gz = os.path.join(tmp_dir, "variants.vcf.gz")
        consensus_pre = os.path.join(tmp_dir, "consensus_pre.fasta")

        # Run mpileup and call variants
        mpileup_cmd = [
            "bcftools", "mpileup",
            "-Ou", "-f", tmp_ref_path,
            bam_path
        ]
        call_cmd = [
            "bcftools", "call",
            "-mv", "-Oz", "-o", vcf_gz
        ]
        try:
            mpileup_process = subprocess.Popen(
                mpileup_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
            call_process = subprocess.run(
                call_cmd, stdin=mpileup_process.stdout,
                check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
            mpileup_process.wait()
            if mpileup_process.returncode != 0:
                raise subprocess.CalledProcessError(
                    mpileup_process.returncode, mpileup_cmd,
                    output=mpileup_process.stdout, stderr=mpileup_process.stderr
                )
        except subprocess.CalledProcessError as e:
            logger.error(f"Variant calling failed: {e.stderr}")
            raise

        # Index VCF
        index_cmd = ["bcftools", "index", vcf_gz]
        try:
            subprocess.run(
                index_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
        except subprocess.CalledProcessError as e:
            logger.error(f"VCF indexing failed: {e.stderr}")
            raise

        # Generate initial consensus
        consensus_cmd = [
            "bcftools", "consensus",
            "-f", tmp_ref_path,
            "-o", consensus_pre,
            vcf_gz
        ]
        try:
            subprocess.run(
                consensus_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
        except subprocess.CalledProcessError as e:
            logger.error(f"Initial consensus failed: {e.stderr}")
            raise

        # Generate low coverage mask BED using bedtools
        bedgraph_all = os.path.join(tmp_dir, "coverage.bedgraph")
        bedgraph_low = os.path.join(tmp_dir, "low_coverage.bedgraph")
        mask_bed = os.path.join(tmp_dir, "mask.bed")

        # Generate full coverage bedgraph
        genomecov_cmd = [
            "bedtools", "genomecov",
            "-ibam", bam_path,
            "-bga"
        ]
        try:
            with open(bedgraph_all, 'w') as f:
                subprocess.run(
                    genomecov_cmd, check=True, stdout=f,
                    stderr=subprocess.PIPE, text=True
                )
        except subprocess.CalledProcessError as e:
            logger.error(f"Genome coverage failed: {e.stderr}")
            raise

        # Filter low coverage regions
        awk_cmd = [
            "awk",
            "-v", f"cov_threshold={cov_thresh}",
            '$4 < cov_threshold {print $1"\t"$2"\t"$3"\t"$4}',
            bedgraph_all
        ]
        try:
            with open(bedgraph_low, 'w') as f:
                subprocess.run(
                    awk_cmd, check=True, stdout=f,
                    stderr=subprocess.PIPE, text=True
                )
        except subprocess.CalledProcessError as e:
            logger.error(f"Low coverage filtering failed: {e.stderr}")
            raise

        # Process and merge BED intervals
        cut_cmd = ["cut", "-f1-3", bedgraph_low]
        sort_cmd = ["bedtools", "sort", "-i", "-"]
        merge_cmd = ["bedtools", "merge", "-i", "-"]

        try:
            cut_process = subprocess.Popen(
                cut_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
            sort_process = subprocess.Popen(
                sort_cmd, stdin=cut_process.stdout,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
            with open(mask_bed, 'w') as f:
                merge_process = subprocess.run(
                    merge_cmd, stdin=sort_process.stdout,
                    check=True, stdout=f, stderr=subprocess.PIPE, text=True
                )
            cut_process.wait()
            sort_process.wait()
            if cut_process.returncode != 0 or sort_process.returncode != 0:
                raise subprocess.CalledProcessError(
                    1, "BED processing pipeline",
                    stderr="BED cutting or sorting failed"
                )
        except subprocess.CalledProcessError as e:
            logger.error(f"BED interval processing failed: {e.stderr}")
            raise

        # Apply mask to initial consensus
        mask_cmd = [
            "bedtools", "maskfasta",
            "-fi", consensus_pre,
            "-bed", mask_bed,
            "-fo", output_fasta
        ]
        try:
            subprocess.run(
                mask_cmd, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True
            )
            logger.info(f"Final masked consensus saved to: {output_fasta}")
        except subprocess.CalledProcessError as e:
            logger.error(f"Consensus masking failed: {e.stderr}")
            raise

    # Clean up temporary reference
    os.unlink(tmp_ref_path)


def main():
    parser = argparse.ArgumentParser(description='Analyze serotype reference alignments and generate consensus sequences')

    required_args = parser.add_argument_group('REQUIRED ARGUMENTS')
    required_args.add_argument("-b", "--bam", required=True, help='Path to sorted BAM file ({sample}_influenza_A.sorted.bam)')
    required_args.add_argument("-s", "--sample", required=True, help='Sample name')
    required_args.add_argument("-o", "--output_dir", required=True, help='Output directory')
    required_args.add_argument("--db", required=True, help='Path to database directory containing Influenza_A_segment_info1.tsv and reference fasta')

    optional_args = parser.add_argument_group('OPTIONAL ARGUMENTS')
    optional_args.add_argument("--cov-thresh", type=int, default=10,
                             help=f"Minimum coverage threshold for masking (default: 10)")

    args = parser.parse_args()

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(os.path.join(args.output_dir, f"{args.sample}_ref_analysis.log")),
            logging.StreamHandler()
        ]
    )

    Path(args.output_dir).mkdir(parents=True, exist_ok=True)

    # Define file paths
    filtered_bam = os.path.join(args.output_dir, f"{args.sample}_filtered.bam")
    coverage_tsv = os.path.join(args.output_dir, f"{args.sample}_coverage.tsv")
    filtered_tsv = os.path.join(args.output_dir, f"{args.sample}_coverage_filtered.tsv")
    annotated_tsv = os.path.join(args.output_dir, f"{args.sample}_coverage_annotated.tsv")
    db_info_path = os.path.join(args.db, "Influenza_A_segment_info1.tsv")
    ref_fasta = os.path.join(args.db, "Influenza_A_segment_sequences.fna")
    vcf_path = os.path.join(args.output_dir, f"{args.sample}_variants.vcf")
    consensus_fasta = os.path.join(args.output_dir, f"{args.sample}_consensus.fasta")

    # Validate database files
    if not os.path.exists(db_info_path):
        logger.error(f"Database info file not found: {db_info_path}")
        return
    if not os.path.exists(ref_fasta):
        logger.error(f"Reference fasta not found: {ref_fasta}")
        return

    # Run analysis pipeline with iav_serotype-compatible filtering
    try:
        logger.info("Filtering BAM with iav_serotype quality criteria...")
        filter_bam_by_quality(args.bam, filtered_bam)

        logger.info("Starting coverage analysis...")
        run_samtools_coverage(filtered_bam, coverage_tsv)

        logger.info("Filtering coverage results...")
        filter_coverage(coverage_tsv, filtered_tsv)

        logger.info("Adding reference genome information...")
        add_reference_info(filtered_tsv, db_info_path, annotated_tsv)

        logger.info(f"Generating consensus sequence with coverage threshold {args.cov_thresh}...")
        generate_consensus(
            filtered_bam, ref_fasta, filtered_tsv,
            vcf_path, consensus_fasta,
            cov_thresh=args.cov_thresh
        )

        logger.info("Analysis completed successfully")
    except Exception as e:
        logger.error(f"Analysis failed: {str(e)}")
        raise


if __name__ == "__main__":
    main()
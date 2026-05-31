#!/usr/bin/env python
"""
Filters a VCF file by excluding duplication SVs (SVTYPE=DUP) that overlap with regions
specified in a BED file.
"""
import logging
from pysam import VariantFile

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("filter_bed_cnvs")

def filter_vcf_to_gz(input_vcf: str, bed_file: str, output_vcf_gz: str) -> None:
    logger.info(f"Loading exclude regions from: {bed_file}")
    exclude_regions = {}
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith(('#', 'track', 'browser')) or not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 3:
                continue
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            exclude_regions.setdefault(chrom, []).append((start, end))

    logger.info(f"Filtering VCF: {input_vcf} -> {output_vcf_gz}")
    vcf_in = VariantFile(input_vcf)
    
    # Auto-detect compressed output format based on extension
    mode = 'wz' if output_vcf_gz.endswith('.gz') else 'w'
    vcf_out = VariantFile(output_vcf_gz, mode, header=vcf_in.header)

    total_records = 0
    excluded_count = 0
    
    for record in vcf_in:
        total_records += 1
        is_dup = record.info.get('SVTYPE') == 'DUP'

        if is_dup:
            chrom = record.chrom
            start, end = record.start, record.stop
            overlap = False
            
            if chrom in exclude_regions:
                for b_start, b_end in exclude_regions[chrom]:
                    if start < b_end and end > b_start:
                        overlap = True
                        break

            if overlap:
                excluded_count += 1
                continue

        vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()
    logger.info(f"Filtering complete! Total records: {total_records}, Excluded DUP overlaps: {excluded_count}")

if __name__ == "__main__":
    filter_vcf_to_gz(snakemake.input.vcf, snakemake.input.bed, snakemake.output.vcf)

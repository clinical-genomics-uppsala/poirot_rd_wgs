import pysam
import sys


def filter_vcf_to_gz(input_vcf, bed_file, output_vcf_gz):
    # 1. Load BED regions
    exclude_regions = {}
    with open(bed_file, 'r') as f:
        for line in f:
            if line.startswith(('#', 'track', 'browser')):
                continue
            parts = line.strip().split('\t')
            chrom, start, end = parts[0], int(parts[1]), int(parts[2])
            exclude_regions.setdefault(chrom, []).append((start, end))

    # 2. Open input and compressed output
    vcf_in = pysam.VariantFile(input_vcf)

    # 'wz' mode specifies 'write' and 'compressed'
    vcf_out = pysam.VariantFile(output_vcf_gz, 'wz', header=vcf_in.header)

    for record in vcf_in:
        is_dup = record.info.get('SVTYPE') == 'DUP'

        if is_dup:
            chrom = record.chrom
            # pysam uses 0-based start (record.start) and 0-based end (record.stop)
            # which aligns perfectly with BED coordinates
            start, end = record.start, record.stop

            overlap = False
            if chrom in exclude_regions:
                for b_start, b_end in exclude_regions[chrom]:
                    if start < b_end and end > b_start:
                        overlap = True
                        break

            if overlap:
                continue

        vcf_out.write(record)

    vcf_in.close()
    vcf_out.close()


if __name__ == "__main__":
    filter_vcf_to_gz(snakemake.input.vcf, snakemake.input.bed, snakemake.output.vcf)

#!/usr/bin/env python
"""
Fixes MT genotype fields in VCF that are longer than 2 alleles (haploid/ploid issues) 
by replacing them with a standard heterozygous call (0, 1).
"""
import logging
from pysam import VariantFile

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("fix_mt_gt")

def fix_gt_field(input_vcf: str, output_vcf: str) -> None:
    logger.info(f"Processing VCF: {input_vcf} -> {output_vcf}")
    vcf_in = VariantFile(input_vcf, 'r')
    
    # Auto-detect compressed output format based on extension
    mode = 'wz' if output_vcf.endswith('.gz') else 'w'
    vcf_out = VariantFile(output_vcf, mode, header=vcf_in.header)
    
    fixed_count = 0
    for variant in vcf_in:
        for sample in variant.samples:
            if 'GT' in variant.samples[sample]:
                gt = variant.samples[sample]['GT']
                if gt and len(gt) > 2:
                    logger.info(f"Fixing GT for sample {sample} at {variant.chrom}:{variant.pos}, {gt} -> (0, 1)")
                    variant.samples[sample]['GT'] = (0, 1)
                    fixed_count += 1
        vcf_out.write(variant)

    vcf_in.close()
    vcf_out.close()
    logger.info(f"Completed! Total genotypes fixed: {fixed_count}")

if __name__ == "__main__":
    fix_gt_field(snakemake.input.vcf, snakemake.output.vcf)

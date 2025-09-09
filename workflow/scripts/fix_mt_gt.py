import pysam


def fix_gt_field(input_vcf, output_vcf):
    
    vcf_in = pysam.VariantFile(input_vcf, 'r')
    
    vcf_out = pysam.VariantFile(output_vcf, 'w', header=vcf_in.header)
    
    for variant in vcf_in:
        for sample in variant.samples:
            # Check if GT field exists and its length
            if 'GT' in variant.samples[sample]:
                gt = variant.samples[sample]['GT']
                if gt and len(gt) > 2:
                    print(f"Fixing GT for sample {sample} at {variant.chrom}:{variant.pos}, {gt} -> (0, 1)")
                    # Replace GT with 0/1
                    variant.samples[sample]['GT'] = (0, 1)
        
        vcf_out.write(variant)
    
    # Close files
    vcf_in.close()
    vcf_out.close()

if __name__ == "__main__":
    
    sample = snakemake.wildcards.sample
    input_vcf = snakemake.input.vcf
    output_vcf = snakemake.output.vcf
    fix_gt_field(input_vcf, output_vcf)
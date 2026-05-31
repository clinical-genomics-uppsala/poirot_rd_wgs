#!/usr/bin/env python
from pysam import VariantFile

def main():
    input_vcf = snakemake.input[0]
    ref_file = snakemake.input[1]
    output_vcf = snakemake.output[0]
    
    vcf_in = VariantFile(input_vcf)
    new_header = vcf_in.header
    new_header.add_line(f"##reference={ref_file}")
    
    # Auto-detect compressed output based on filename extension
    mode = 'wz' if output_vcf.endswith('.gz') else 'w'
    
    vcf_out = VariantFile(output_vcf, mode, header=new_header)
    for record in vcf_in.fetch():
        vcf_out.write(record)
        
    vcf_in.close()
    vcf_out.close()

if __name__ == "__main__":
    main()

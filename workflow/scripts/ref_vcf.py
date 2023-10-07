#!/bin/python3.6

import gzip

ref_str = snakemake.input[1]
with gzip.open(snakemake.input[0], 'r') as vcf_in:
    with open(snakemake.output[0], 'w') as vcf_out:
        line_count = 0
        for line in vcf_in:
            line_count += 1
            if line_count == 4:
                print(f"##reference={ref_str}", file=vcf_out)
            print(line.decode(encoding="utf-8").rstrip(), file=vcf_out)

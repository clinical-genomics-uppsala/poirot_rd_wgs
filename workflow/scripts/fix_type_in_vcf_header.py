
with open(snakemake.input.vcf, 'r') as vcf_in:
    with open(snakemake.output[0], 'w') as vcf_out:
        for line in vcf_in:
            if line.startswith("##INFO=<ID=gnomad_AF"):
                col = line.rstrip().split(',')
                col[2] = 'Type=Float'
                print(','.join(col), file=vcf_out)
            else:
                print(line.rstrip(), file=vcf_out)
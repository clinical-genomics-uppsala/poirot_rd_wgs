#!/usr/bin/env python

from pysam import VariantFile

vcf_in = snakemake.input.vcf
aed = snakemake.output.aed

vcf = VariantFile(vcf_in)

header_0 = ["bio:sequence(aed:String)", "bio:start(aed:Integer)",
            "bio:end(aed:Integer)", "aed:name(aed:String)",
            "aed:value(aed:String)", "aed:category(aed:String)",
            "style:color(aed:Color)"]

header_1 = ["", "", "", "namespace:affx(aed:URI)",
            "http://affymetrix.com/ontology/", "", ""]

header_2 = ["", "", "", "affx:ucscGenomeVersion(aed:String)", "hg38", "", ""]


with open(aed, 'w') as outfile:
    print('\t'.join(header_0), file=outfile)
    print('\t'.join(header_1), file=outfile)
    print('\t'.join(header_2), file=outfile)
    for rec in vcf:
        var_type = rec.info["SVTYPE"]
        if var_type == 'DEL':
            aed_vartype = "copynumber/loss"
            colour = "rgb(255,0,0)"
        elif var_type == 'DUP':
            aed_vartype = "copynumber/gain"
            colour = "rgb(0,0,255)"

        print(rec.contig, rec.start,
              rec.stop, rec.id, "", aed_vartype,
              colour, sep="\t", file=outfile)

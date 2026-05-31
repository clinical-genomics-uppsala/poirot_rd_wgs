#!/usr/bin/env python
"""
Converts CNVPytor VCF output into AED (Affymetrix Event Detection) format.
"""
import logging
from pysam import VariantFile

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("cnvpytor_vcf_to_aed")

def convert_vcf_to_aed(input_vcf: str, output_aed: str) -> None:
    logger.info(f"Converting CNVPytor VCF to AED: {input_vcf} -> {output_aed}")
    vcf = VariantFile(input_vcf)

    header_0 = ["bio:sequence(aed:String)", "bio:start(aed:Integer)",
                "bio:end(aed:Integer)", "aed:name(aed:String)",
                "aed:value(aed:String)", "aed:category(aed:String)",
                "style:color(aed:Color)"]

    header_1 = ["", "", "", "namespace:affx(aed:URI)",
                "http://affymetrix.com/ontology/", "", ""]

    header_2 = ["", "", "", "affx:ucscGenomeVersion(aed:String)", "hg38", "", ""]

    record_count = 0
    with open(output_aed, 'w') as outfile:
        print('\t'.join(header_0), file=outfile)
        print('\t'.join(header_1), file=outfile)
        print('\t'.join(header_2), file=outfile)
        
        for rec in vcf:
            var_type = rec.info.get("SVTYPE")
            if var_type == 'DEL':
                aed_vartype = "copynumber/loss"
                colour = "rgb(255,0,0)"
            elif var_type == 'DUP':
                aed_vartype = "copynumber/gain"
                colour = "rgb(0,0,255)"
            else:
                continue

            print(rec.contig, rec.start, rec.stop, rec.id or "", "", 
                  aed_vartype, colour, sep="\t", file=outfile)
            record_count += 1

    vcf.close()
    logger.info(f"Conversion complete! Total events exported: {record_count}")

if __name__ == "__main__":
    convert_vcf_to_aed(snakemake.input.vcf, snakemake.output.aed)

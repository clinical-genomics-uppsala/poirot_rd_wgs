import pysam


loci_to_filter = ["HTT", "HTT_CCG", "C9ORF72"]
vcf = pysam.VariantFile(snakemake.input.vcf)
if snakemake.input.panel_list is not None:
    loci_to_keep = [locus.rstrip() for locus in open(
        snakemake.input.panel_list, 'r')]


def get_bed_rec(rec):
    try:
        str_status = rec.info["STR_STATUS"]
        str_normal_max = str(rec.info["STR_NORMAL_MAX"])
        str_pathologic_min = str(rec.info["STR_PATHOLOGIC_MIN"])
    except KeyError:
        str_status = "NA"
        str_normal_max = "NA"
        str_pathologic_min = "NA"

    if rec.alts is not None:
        alts = ",".join(rec.alts)
    else:
        alts = "."

    if str_status != "NA":
        str_status = ','.join(str_status)

    bed_rec = "\t".join([rec.chrom, str(rec.pos - 1), str(rec.stop),
                        rec.info["VARID"], alts,
                        str_status, str_normal_max,
                        str_pathologic_min])

    return bed_rec


with open(snakemake.output.bed, "w") as outfile:
    print("#chrom\tstart\tend\tvariant_id\talt\t\
          str_status\tstr_normal_max\t\
          str_pathologic_min", file=outfile)
    for rec in vcf:
        bed_rec = get_bed_rec(rec)
        if snakemake.input.panel_list is not None:
            if rec.info["VARID"] in loci_to_keep:
                print(bed_rec, file=outfile)
        else:
            if rec.info["VARID"] not in loci_to_filter:
                print(bed_rec, file=outfile)

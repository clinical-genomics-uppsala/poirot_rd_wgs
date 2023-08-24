#!/usr/bin/env python

vcf_in = snakemake.input.vcf
aed = snakemake.output.aed

def get_info_dict(info_list):
    info_dict = {}
    for i in info_list:
        info_rec = i.split("=")
        if info_rec[0] == 'IMPRECISE':
            continue
        info_dict[info_rec[0]] = info_rec[1]

    return info_dict


header_0 = ["bio:sequence(aed:String)", "bio:start(aed:Integer)",
            "bio:end(aed:Integer)", "aed:name(aed:String)",
            "aed:value(aed:String)","aed:category(aed:String)",
            "style:color(aed:Color)"]

header_1 = ["", "", "", "namespace:affx(aed:URI)",
            "http://affymetrix.com/ontology/", "", ""]

header_2 = ["", "", "", "affx:ucscGenomeVersion(aed:String)", "hg38", "", ""]


with open(aed, 'w') as outfile:
    print('\t'.join(header_0), file=outfile)
    print('\t'.join(header_1), file=outfile)
    print('\t'.join(header_2), file=outfile)
    with open(vcf_in, 'r') as vcf:
        for rec in vcf:
            if rec.startswith('#'):
                continue
            vcf_rec = rec.split('\t')
            contig = vcf_rec[0]
            start = int(vcf_rec[1]) - 1
            rec_id = vcf_rec[2]
            info_list = vcf_rec[7].split(";")
            info_dict = get_info_dict(info_list)
            var_type = info_dict["SVTYPE"]
            end = info_dict["END"]
            if var_type == 'DEL':
                aed_vartype = "copynumber/loss"
                colour = "rgb(255,0,0)"
            elif var_type == 'DUP':
                aed_vartype = "copynumber/gain"
                colour = "rgb(0,0,255)"

            print(contig, start, end, rec_id, "", aed_vartype, colour,
                  sep="\t", file=outfile)

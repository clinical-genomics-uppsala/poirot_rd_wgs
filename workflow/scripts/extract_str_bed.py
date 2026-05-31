#!/usr/bin/env python
"""
Extracts STR loci from VCF and writes them as a BED file, 
optionally filtering using a whitelist or blacklist.
"""
import logging
import pysam

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
logger = logging.getLogger("extract_str_bed")

def extract_str(vcf_file: str, bed_output: str, panel_list_file: str = None) -> None:
    loci_to_filter = {"HTT", "HTT_CCG", "C9ORF72"}
    
    loci_to_keep = []
    if panel_list_file:
        try:
            with open(panel_list_file, 'r') as f:
                for line in f:
                    if line.strip():
                        loci_to_keep.append(line.strip())
            logger.info(f"Loaded {len(loci_to_keep)} loci to keep from whitelist: {panel_list_file}")
        except Exception as e:
            logger.warning(f"Could not read panel list whitelist '{panel_list_file}': {e}")

    logger.info(f"Extracting STR loci: {vcf_file} -> {bed_output}")
    vcf = pysam.VariantFile(vcf_file)

    def get_bed_rec(rec):
        try:
            str_status = rec.info.get("STR_STATUS", "NA")
            str_normal_max = str(rec.info.get("STR_NORMAL_MAX", "NA"))
            str_pathologic_min = str(rec.info.get("STR_PATHOLOGIC_MIN", "NA"))
        except KeyError:
            str_status = "NA"
            str_normal_max = "NA"
            str_pathologic_min = "NA"

        alts = ",".join(rec.alts) if rec.alts else "."
        
        if str_status != "NA" and not isinstance(str_status, str):
            str_status = ','.join(str_status)

        var_id = rec.info.get("VARID", "NA")
        
        # pysam uses 0-based start, which perfectly matches BED coordinates
        return "\t".join([
            rec.chrom, 
            str(rec.pos - 1), 
            str(rec.stop),
            var_id, 
            alts,
            str_status, 
            str_normal_max,
            str_pathologic_min
        ])

    written_count = 0
    with open(bed_output, "w") as outfile:
        header = "#chrom\tstart\tend\tvariant_id\talt\tstr_status\tstr_normal_max\tstr_pathologic_min"
        print(header, file=outfile)
        
        for rec in vcf:
            var_id = rec.info.get("VARID")
            if not var_id:
                continue
                
            bed_rec = get_bed_rec(rec)
            if len(loci_to_keep) > 0:
                if var_id in loci_to_keep:
                    print(bed_rec, file=outfile)
                    written_count += 1
            else:
                if var_id not in loci_to_filter:
                    print(bed_rec, file=outfile)
                    written_count += 1

    vcf.close()
    logger.info(f"STR extraction complete! Total records written: {written_count}")

if __name__ == "__main__":
    # Check if panel_list is provided as an input attribute
    panel_list = None
    try:
        panel_list = snakemake.input.panel_list
    except AttributeError:
        pass
        
    extract_str(snakemake.input.vcf, snakemake.output.bed, panel_list)

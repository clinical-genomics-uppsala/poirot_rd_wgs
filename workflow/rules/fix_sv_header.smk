__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2024, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule fix_sv_header:
    input:
        vcf="annotate/vep_svdb/{sample}_{type}.merged.svdb_query.vep_annotated.vcf",
    output:
        vcf="annotate/vep_svdb/{sample}_{type}.merged.svdb_query.vep_annotated.fixed_header.vcf",
    params:
        extra=config.get("fix_sv_header", {}).get("extra", ""),
    log:
        "poirot_rd_wgs/fix_sv_header/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "poirot_rd_wgs/fix_sv_header/{sample}_{type}.output.benchmark.tsv",
            config.get("fix_sv_header", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("fix_sv_header", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fix_sv_header", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fix_sv_header", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("fix_sv_header", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("fix_sv_header", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fix_sv_header", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fix_sv_header", {}).get("container", config["default_container"])
    message:
        "{rule}: remove the gnomad info headers added by vep --custom due to the type being always string {input.vcf}"
    script:
        "../scripts/fix_type_in_vcf_header.py"

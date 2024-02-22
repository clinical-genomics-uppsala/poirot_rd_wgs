__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2024, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"



rule bcftools_split_vep:
    input:
        vcf="annotate/vep_svdb/{sample}_{type}.merged.svdb_query.vep_annotated.fixed_header.vcf",
    output:
        vcf="annotate/vep_svdb/{sample}_{type}.merged.svdb_query.vep_annotated.vep_info.vcf",
    params:
        columns=config.get("bcftools_split_vep", {}).get("columns", "gnomad_AF:Float"),
        extra=config.get("bcftools_split_vep", {}).get("extra", ""),
    log:
        "poirot_rd_wgs/bcftools_split_vep/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "poirot_rd_wgs/bcftools_split_vep/{sample}_{type}.output.benchmark.tsv",
            config.get("bcftools_split_vep", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("bcftools_split_vep", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bcftools_split_vep", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcftools_split_vep", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcftools_split_vep", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcftools_split_vep", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bcftools_split_vep", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcftools_split_vep", {}).get("container", config["default_container"])
    message:
        "{rule}: add vep annotations to the info field {input.vcf}"
    shell:
        "(bcftools +split-vep "
        "-c {params.columns} "
        "{input.vcf} "
        " > {output.vcf}) &> {log}"

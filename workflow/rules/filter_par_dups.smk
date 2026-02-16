__author__ = "Pádraic Corcoran"
__copyright__ = "Copyright 2026, Pádraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule filter_par_dups:
    input:
        vcf=get_cnvpytor_male_input,
        bed=config["filter_par_dups"]["bed"],
    output:
        vcf="cnv_sv/cnvpytor/{sample}_{type}.par_dups_filtered.vcf.gz",
    params:
        extra=config.get("filter_par_dups", {}).get("extra", ""),
    log:
        "cnv_sv/cnvpytor/{sample}_{type}.par_dups_filtered.vcf.gz.log",
    benchmark:
        repeat(
            "cnv_sv/cnvpytor/{sample}_{type}.par_dups_filtered.vcf.gz.benchmark.tsv",
            config.get("filter_par_dups", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("filter_par_dups", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("filter_par_dups", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("filter_par_dups", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("filter_par_dups", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("filter_par_dups", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("filter_par_dups", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("filter_par_dups", {}).get("container", config["default_container"])
    message:
        "{rule}: filter cnvpytor DUP calls in {input.vcf} located in for {input.bed}"
    script:
        "../scripts/filter_bed_cnvs.py"

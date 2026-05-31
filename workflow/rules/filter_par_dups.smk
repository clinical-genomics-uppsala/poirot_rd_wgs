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
    threads: rule_resource("filter_par_dups", "threads")
    resources:
        mem_mb=rule_resource("filter_par_dups", "mem_mb"),
        mem_per_cpu=rule_resource("filter_par_dups", "mem_per_cpu"),
        partition=rule_resource("filter_par_dups", "partition"),
        threads=rule_resource("filter_par_dups", "threads"),
        time=rule_resource("filter_par_dups", "time"),
    container:
        rule_container("filter_par_dups")
    message:
        "{rule}: filter cnvpytor DUP calls in {input.vcf} located in for {input.bed}"
    script:
        "../scripts/filter_bed_cnvs.py"

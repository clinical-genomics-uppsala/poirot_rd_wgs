__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2023, Padraic Corcoran"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule vcf_to_aed:
    input:
        vcf="cnv_sv/cnvpytor/{sample}_{type}.vcf",
    output:
        aed="cnv_sv/cnvpytor/{sample}_{type}.aed",
    params:
        extra=config.get("vcf_to_aed", {}).get("extra", ""),
    log:
        "cnv_sv/cnvpytor/{sample}_{type}.aed.log",
    benchmark:
        repeat("cnv_sv/cnvpytor/{sample}_{type}.aed.benchmark.tsv", config.get("vcf_to_aed", {}).get("benchmark_repeats", 1))
    threads: rule_resource("vcf_to_aed", "threads")
    resources:
        mem_mb=rule_resource("vcf_to_aed", "mem_mb"),
        mem_per_cpu=rule_resource("vcf_to_aed", "mem_per_cpu"),
        partition=rule_resource("vcf_to_aed", "partition"),
        threads=rule_resource("vcf_to_aed", "threads"),
        time=rule_resource("vcf_to_aed", "time"),
    container:
        rule_container("vcf_to_aed")
    message:
        "{rule}: convert {input.vcf} to AED format"
    script:
        "../scripts/cnvpytor_vcf_to_aed.py"


rule vcf_to_aed_filtered:
    input:
        vcf="cnv_sv/cnvpytor/{sample}_{type}.hardfiltered.vcf",
    output:
        aed="cnv_sv/cnvpytor/{sample}_{type}_filtered.aed",
    params:
        extra=config.get("vcf_to_aed", {}).get("extra", ""),
    log:
        "cnv_sv/cnvpytor/{sample}_{type}.aed.log",
    benchmark:
        repeat("cnv_sv/cnvpytor/{sample}_{type}.aed.benchmark.tsv", config.get("vcf_to_aed", {}).get("benchmark_repeats", 1))
    threads: rule_resource("vcf_to_aed", "threads")
    resources:
        mem_mb=rule_resource("vcf_to_aed", "mem_mb"),
        mem_per_cpu=rule_resource("vcf_to_aed", "mem_per_cpu"),
        partition=rule_resource("vcf_to_aed", "partition"),
        threads=rule_resource("vcf_to_aed", "threads"),
        time=rule_resource("vcf_to_aed", "time"),
    container:
        rule_container("vcf_to_aed")
    message:
        "{rule}: convert {input.vcf} to AED format"
    script:
        "../scripts/cnvpytor_vcf_to_aed.py"

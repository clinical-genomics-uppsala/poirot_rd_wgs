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
    threads: config.get("vcf_to_aed", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("vcf_to_aed", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("vcf_to_aed", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("vcf_to_aed", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("vcf_to_aed", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("vcf_to_aed", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("vcf_to_aed", {}).get("container", config["default_container"])
    message:
        "{rule}: convert {input.vcf} to AED format"
    script:
        "../scripts/cnvpytor_vcf_to_aed.py"


rule vcf_to_aed_filtered:
    input:
        vcf="cnv_sv/cnvpytor/{sample}_{type}.filtered.vcf",
    output:
        aed="cnv_sv/cnvpytor/{sample}_{type}_filtered.aed",
    params:
        extra=config.get("vcf_to_aed", {}).get("extra", ""),
    log:
        "cnv_sv/cnvpytor/{sample}_{type}.aed.log",
    benchmark:
        repeat("cnv_sv/cnvpytor/{sample}_{type}.aed.benchmark.tsv", config.get("vcf_to_aed", {}).get("benchmark_repeats", 1))
    threads: config.get("vcf_to_aed", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("vcf_to_aed", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("vcf_to_aed", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("vcf_to_aed", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("vcf_to_aed", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("vcf_to_aed", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("vcf_to_aed", {}).get("container", config["default_container"])
    message:
        "{rule}: convert {input.vcf} to AED format"
    script:
        "../scripts/cnvpytor_vcf_to_aed.py"

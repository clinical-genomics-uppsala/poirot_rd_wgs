__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2022"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule exclude_chrM:
    input:
        deepvariant_vcf="vcf_final/{sample}_ref.vcf",
    output:
        vcf=temp("vcf_final/{sample}.no_ChrM.vcf"),
    log:
        "vcf_final/{sample}.no_ChrM.vcf.log",
    benchmark:
        repeat(
            "vcf_final/{sample}.no_ChrM.vcf.benchmark.tsv",
            config.get("exclude_chrM", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("exclude_chrM", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("exclude_chrM", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("exclude_chrM", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("exclude_chrM", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("exclude_chrM", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("exclude_chrM", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("exclude_chrM", {}).get("container", config["default_container"])
    conda:
        "../envs/concat_vcfs.yaml"
    message:
        "{rule}: Exclude chrM from the deepvariant vcf: {input.deepvariant_vcf}"
    shell:
        "(bcftools view "
        "-t ^chrM "
        "{input.deepvariant_vcf} "
        "-o {output.vcf}) &> {log}"


rule bcftools_concat:
    input:
        deepvariant_vcf="vcf_final/{sample}.no_ChrM.vcf",
        mutect2_vcf="mitochondrial/gatk_select_variants_final/{sample}_N.vcf",
    output:
        vcf=temp("vcf_final/{sample}.vcf"),
    log:
        "vcf_final/{sample}.vcf.log",
    benchmark:
        repeat(
            "vcf_final/{sample}.vcf.benchmark.tsv",
            config.get("bcftools_concat", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("bcftools_concat", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("bcftools_concat", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcftools_concat", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcftools_concat", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcftools_concat", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bcftools_concat", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcftools_concat", {}).get("container", config["default_container"])
    conda:
        "../envs/concat_vcfs.yaml"
    message:
        "{rule}: Concatenate the deepvariant vcf: {input.deepvariant_vcf} and mitochondrial vcf: {input.mutect2_vcf}"
    shell:
        "(bcftools concat "
        "{input.deepvariant_vcf} "
        "{input.mutect2_vcf} "
        "-o {output.vcf}) &> {log}"

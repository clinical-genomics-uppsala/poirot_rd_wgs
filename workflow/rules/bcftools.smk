__author__ = "Jessika Nordin"
__copyright__ = "Copyright 2023"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule bcftools_pass:
    input:
        vcf="vcf_final/{sample}_{type}.vcf.gz",
    output:
        vcf=temp("vcf_final/{sample}_{type}.pass.vcf.gz"),
    log:
        "vcf_final/{sample}_{type}.pass.vcf.log",
    benchmark:
        repeat(
            "vcf_final/{sample}_{type}.pass.vcf.benchmark.tsv",
            config.get("pass", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("pass", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("pass", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("pass", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("pass", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("pass", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("pass", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("pass", {}).get("container", config["default_container"])
    conda:
        "../envs/bcftools.yaml"
    message:
        "{rule}: Only keep variants with PASS in FILTER column from the deepvariant vcf: {input.vcf}"
    shell:
        "(bcftools view "
        "--apply-filters PASS "
        "{input.vcf} "
        "--output-type z --output-file {output.vcf}) &> {log}"

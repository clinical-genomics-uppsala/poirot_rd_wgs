__author__ = "Pádraic Corcoran"
__copyright__ = "Copyright 2025, Pádraic Corcoran"
__email__ = "padraicc"
__license__ = "GPL-3"


rule fix_mt_gt:
    input:
        vcf="mitochondrial/gatk_select_variants_final/{sample}_{type}.fix_af.vcf",
    output:
        vcf="mitochondrial/gatk_select_variants_final/{sample}_{type}.fix_gt.vcf",
    params:
        extra=config.get("fix_mt_gt", {}).get("extra", ""),
    log:
        "mitochondrial/gatk_select_variants_final/{sample}_{type}.fix_af.gt_fixed.vcf",
    benchmark:
        repeat(
            "mitochondrial/gatk_select_variants_final/{sample}_{type}.fix_af.gt_fixed.vcf.benchmark.tsv",
            config.get("fix_mt_gt", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("fix_mt_gt", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("fix_mt_gt", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("fix_mt_gt", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("fix_mt_gt", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("fix_mt_gt", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("fix_mt_gt", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("fix_mt_gt", {}).get("container", config["default_container"])
    message:
        "{rule}: fix GT fields with >2 alleles in {input.vcf} after GATK multiallelic splitting (e.g., '0/././1')."
    script:
        "../scripts/fix_mt_gt.py"

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
    threads: rule_resource("fix_mt_gt", "threads")
    resources:
        mem_mb=rule_resource("fix_mt_gt", "mem_mb"),
        mem_per_cpu=rule_resource("fix_mt_gt", "mem_per_cpu"),
        partition=rule_resource("fix_mt_gt", "partition"),
        threads=rule_resource("fix_mt_gt", "threads"),
        time=rule_resource("fix_mt_gt", "time"),
    container:
        rule_container("fix_mt_gt")
    message:
        "{rule}: fix GT fields with >2 alleles in {input.vcf} after GATK multiallelic splitting (e.g., '0/././1')."
    script:
        "../scripts/fix_mt_gt.py"

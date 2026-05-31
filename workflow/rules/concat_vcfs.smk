__author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2022"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule exclude_chrM:
    input:
        deepvariant_vcf="snv_indels/vcf_final/{sample}_{type}_ref.vcf",
    output:
        vcf=temp("snv_indels/vcf_final/{sample}_{type}.no_ChrM.vcf"),
    log:
        "snv_indels/vcf_final/{sample}_{type}.no_ChrM.vcf.log",
    benchmark:
        repeat(
            "snv_indels/vcf_final/{sample}_{type}.no_ChrM.vcf.benchmark.tsv",
            config.get("exclude_chrM", {}).get("benchmark_repeats", 1),
        )
    threads: rule_resource("exclude_chrM", "threads")
    resources:
        mem_mb=rule_resource("exclude_chrM", "mem_mb"),
        mem_per_cpu=rule_resource("exclude_chrM", "mem_per_cpu"),
        partition=rule_resource("exclude_chrM", "partition"),
        threads=rule_resource("exclude_chrM", "threads"),
        time=rule_resource("exclude_chrM", "time"),
    container:
        rule_container("exclude_chrM")
    message:
        "{rule}: Exclude chrM from the deepvariant vcf: {input.deepvariant_vcf}"
    shell:
        "(bcftools view "
        "-t ^chrM "
        "{input.deepvariant_vcf} "
        "-o {output.vcf}) &> {log}"


rule bcftools_concat:
    input:
        deepvariant_vcf="snv_indels/vcf_final/{sample}_{type}.no_ChrM.vcf",
        mutect2_vcf="mitochondrial/gatk_select_variants_final/{sample}_{type}.fix_gt.vcf",
    output:
        vcf=temp("snv_indels/vcf_final/{sample}_{type}.vcf"),
    log:
        "snv_indels/vcf_final/{sample}_{type}.vcf.log",
    benchmark:
        repeat(
            "snv_indels/vcf_final/{sample}_{type}.vcf.benchmark.tsv",
            config.get("bcftools_concat", {}).get("benchmark_repeats", 1),
        )
    threads: rule_resource("bcftools_concat", "threads")
    resources:
        mem_mb=rule_resource("bcftools_concat", "mem_mb"),
        mem_per_cpu=rule_resource("bcftools_concat", "mem_per_cpu"),
        partition=rule_resource("bcftools_concat", "partition"),
        threads=rule_resource("bcftools_concat", "threads"),
        time=rule_resource("bcftools_concat", "time"),
    container:
        rule_container("bcftools_concat")
    message:
        "{rule}: Concatenate the deepvariant vcf: {input.deepvariant_vcf} and mitochondrial vcf: {input.mutect2_vcf}"
    shell:
        "(bcftools concat "
        "{input.deepvariant_vcf} "
        "{input.mutect2_vcf} "
        "-o {output.vcf}) &> {log}"

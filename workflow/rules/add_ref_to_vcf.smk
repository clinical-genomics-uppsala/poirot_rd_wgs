__author__ = "Jessika Nordin, Padraic Corcoran"
__copyright__ = "Copyright 2022"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule deepvariant_add_ref:
    input:
        vcf="snv_indels/vcf_final/{sample}_{type}.fix_af.vcf.gz",
        ref=config["reference"]["fasta"],
    output:
        vcf=temp("snv_indels/vcf_final/{sample}_{type}_ref.vcf"),
    log:
        "snv_indels/vcf_final/{sample}_{type}_ref.log",
    benchmark:
        repeat(
            "snv_indels/vcf_final/{sample}_{type}_ref.vcf.benchmark.tsv",
            config.get("deepvariant_add_ref", {}).get("benchmark_repeats", 1),
        )
    resources:
        mem_mb=rule_resource("deepvariant_add_ref", "mem_mb"),
        mem_per_cpu=rule_resource("deepvariant_add_ref", "mem_per_cpu"),
        partition=rule_resource("deepvariant_add_ref", "partition"),
        threads=rule_resource("deepvariant_add_ref", "threads"),
        time=rule_resource("deepvariant_add_ref", "time"),
    container:
        rule_container("deepvariant_add_ref")
    message:
        "{rule}: Add reference to the header of the deepvariant vcf: {input.vcf}"
    script:
        "../scripts/ref_vcf.py"


rule tiddit_add_ref:
    input:
        vcf="cnv_sv/tiddit/{sample}_{type}.vcf.gz",
        ref=config["reference"]["fasta"],
    output:
        vcf="cnv_sv/tiddit/{sample}_{type}_ref.vcf",
    params:
        extra=config.get("tiddit_add_ref", {}).get("extra", ""),
    log:
        "cnv_sv/tiddit/{sample}_{type}.add_ref.log",
    benchmark:
        repeat(
            "cnv_sv/tiddit/{sample}_{type}.add_ref.benchmark.tsv",
            config.get("tiddit_add_ref", {}).get("benchmark_repeats", 1),
        )
    threads: rule_resource("tiddit_add_ref", "threads")
    resources:
        mem_mb=rule_resource("tiddit_add_ref", "mem_mb"),
        mem_per_cpu=rule_resource("tiddit_add_ref", "mem_per_cpu"),
        partition=rule_resource("tiddit_add_ref", "partition"),
        threads=rule_resource("tiddit_add_ref", "threads"),
        time=rule_resource("tiddit_add_ref", "time"),
    container:
        rule_container("tiddit_add_ref")
    message:
        "{rule}: Add reference to the header of the tiddit vcf: {input.vcf}"
    script:
        "../scripts/ref_vcf.py"


rule svdb_add_ref:
    input:
        vcf="cnv_sv/svdb_query/{sample}_{type}.merged.svdb_query.vcf.gz",
        ref=config["reference"]["fasta"],
    output:
        vcf="cnv_sv/svdb_query/{sample}_{type}.merged.svdb_query_ref.vcf",
    params:
        extra=config.get("svdb_add_ref", {}).get("extra", ""),
    log:
        "cnv_sv/svdb_query/{sample}_{type}.add_ref.log",
    benchmark:
        repeat(
            "cnv_sv/svdb_merge/{sample}_{type}.output.benchmark.tsv", config.get("svdb_add_ref", {}).get("benchmark_repeats", 1)
        )
    threads: rule_resource("svdb_add_ref", "threads")
    resources:
        mem_mb=rule_resource("svdb_add_ref", "mem_mb"),
        mem_per_cpu=rule_resource("svdb_add_ref", "mem_per_cpu"),
        partition=rule_resource("svdb_add_ref", "partition"),
        threads=rule_resource("svdb_add_ref", "threads"),
        time=rule_resource("svdb_add_ref", "time"),
    container:
        rule_container("svdb_add_ref")
    message:
        "{rule}: Add reference to the header of the svdb vcf: {input.vcf}"
    script:
        "../scripts/ref_vcf.py"

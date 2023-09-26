__author__ = "Jessika Nordin, Padraic Corcoran"
__copyright__ = "Copyright 2022"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule deepvariant_add_ref:
    input:
        vcf="vcf_final/{sample}_{type}.fix_af.vcf",
        ref=config["reference"]["fasta"],
    output:
        vcf=temp("vcf_final/{sample}_{type}_ref.vcf"),
    log:
        "vcf_final/{sample}_{type}_ref.log",
    benchmark:
        repeat(
            "vcf_final/{sample}_{type}_ref.vcf.benchmark.tsv",
            config.get("deepvariant_add_ref", {}).get("benchmark_repeats", 1),
        )
    resources:
        mem_mb=config.get("deepvariant_add_ref", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("deepvariant_add_ref", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("deepvariant_add_ref", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("deepvariant_add_ref", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("deepvariant_add_ref", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("deepvariant_add_ref", {}).get("container", config["default_container"])
    conda:
        "../envs/deepvariant_add_ref.yaml"
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
        "poirot_rd_wgs/tiddit_add_ref/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "poirot_rd_wgs/tiddit_add_ref/{sample}_{type}.output.benchmark.tsv",
            config.get("tiddit_add_ref", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("tiddit_add_ref", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("tiddit_add_ref", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("tiddit_add_ref", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("tiddit_add_ref", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("tiddit_add_ref", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("tiddit_add_ref", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("tiddit_add_ref", {}).get("container", config["default_container"])
    message:
        "{rule}: Add reference to the header of the tiddit vcf: {input.vcf}"
    script:
        "../scripts/ref_vcf.py"


rule svdb_add_ref:
    input:
        vcf="cnv_sv/svdb_merge/{sample}_{type}.merged.vcf.gz",
        ref=config["reference"]["fasta"],
    output:
        vcf="cnv_sv/svdb_merge/{sample}_{type}.merged_ref.vcf",
    params:
        extra=config.get("svdb_add_ref", {}).get("extra", ""),
    log:
        "cnv_sv/svdb_merge/{sample}_{type}.output.log",
    benchmark:
        repeat(
            "cnv_sv/svdb_merge/{sample}_{type}.output.benchmark.tsv",
            config.get("svdb_add_ref", {}).get("benchmark_repeats", 1)
        )
    threads: config.get("svdb_add_ref", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("svdb_add_ref", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("svdb_add_ref", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("svdb_add_ref", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("svdb_add_ref", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("svdb_add_ref", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("svdb_add_ref", {}).get("container", config["default_container"])
    message:
        "{rule}: Add reference to the header of the svdb vcf: {input.vcf}"
    script:
        "../scripts/ref_vcf.py"

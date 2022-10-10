
rule addRef:
    input:
        vcf="parabricks/pbrun_deepvariant/{sample}.vcf",
        ref=config["reference"]["fasta"],
    output:
        temp("vcf_final/{sample}_ref.vcf"),
    log:
        "vcf_final/{sample}_add_ref.log",
    params:
        config["programdir"]["dir"],
    resources:
        mem_mb=config.get("multiqc", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("multiqc", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("multiqc", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("multiqc", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("multiqc", {}).get("time", config["default_resources"]["time"]),
    shell:
        "( python {params}/scripts/ref_vcf.py {input.vcf} {input.ref} {output} ) &> {log}"


rule changeM2MT:
    input:
        "vcf_final/{sample}_ref.vcf",
    output:
        temp("vcf_final/{sample}.vcf"),
    log:
        "vcf_final/{sample}_chrMT.log",
    resources:
        mem_mb=config.get("multiqc", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("multiqc", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("multiqc", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("multiqc", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("multiqc", {}).get("time", config["default_resources"]["time"]),
    shell:
        """( awk '{{gsub(/chrM/,"chrMT"); print}}' {input} > {output} ) &> {log}"""


# rule bgzipNtabix:
#     input:
#         "vcf_final/{sample}.vcf",
#     output:
#         "vcf_final/{sample}.vcf.gz",
#         "vcf_final/{sample}.vcf.gz.tbi",
#     log:
#         "vcf_final/{sample}.bgzip-tabix.log",
#     resources:
#         mem_mb=config.get("multiqc", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
#         mem_per_cpu=config.get("multiqc", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
#         partition=config.get("multiqc", {}).get("partition", config["default_resources"]["partition"]),
#         threads=config.get("multiqc", {}).get("threads", config["default_resources"]["threads"]),
#         time=config.get("multiqc", {}).get("time", config["default_resources"]["time"]),
#     shell:
#         "( bgzip {input} && tabix {input}.gz ) &> {log}"

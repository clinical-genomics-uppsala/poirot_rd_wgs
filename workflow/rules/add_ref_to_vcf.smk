___author__ = "Jessika Nordin"
__copyright__ = "Copyright 2022"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule vcf_addRef:
    input:
        vcf="parabricks/pbrun_deepvariant/{sample}.vcf",
        ref=config["reference"]["fasta"],
    output:
        vcf=temp("vcf_final/{sample}.vcf"),
    log:
        "vcf_final/{sample}_add_ref.log",
    resources:
        mem_mb=config.get("multiqc", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("multiqc", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("multiqc", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("multiqc", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("multiqc", {}).get("time", config["default_resources"]["time"]),
    script:
        "../scripts/ref_vcf.py"


# vcf=temp("vcf_final/{sample}_ref.vcf"),

# rule vcf_changeM2MT:
#     input:
#         "vcf_final/{sample}_ref.vcf",
#     output:
#         vcf=temp("vcf_final/{sample}.vcf"),
#     log:
#         "vcf_final/{sample}_chrMT.log",
#     resources:
#         mem_mb=config.get("multiqc", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
#         mem_per_cpu=config.get("multiqc", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
#         partition=config.get("multiqc", {}).get("partition", config["default_resources"]["partition"]),
#         threads=config.get("multiqc", {}).get("threads", config["default_resources"]["threads"]),
#         time=config.get("multiqc", {}).get("time", config["default_resources"]["time"]),
#     shell:
#         """( awk '{{gsub(/chrM/,"chrMT"); print}}' {input} > {output.vcf} ) &> {log}"""

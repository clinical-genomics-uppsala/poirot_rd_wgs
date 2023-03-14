__author__ = "Jessika Nordin, Padraic Corcoran"
__copyright__ = "Copyright 2022"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule vcf_addRef:
    input:
        vcf="parabricks/pbrun_deepvariant/{sample}_{type}.fix_af.vcf",
        ref=config["reference"]["fasta"],
    output:
        vcf=temp("vcf_final/{sample}_{type}_ref.vcf"),
    log:
        "vcf_final/{sample}_{type}_ref.log",
    benchmark:
        repeat(
            "vcf_final/{sample}_{type}_ref.vcf.benchmark.tsv",
            config.get("vcf_addRef", {}).get("benchmark_repeats", 1),
        )
    resources:
        mem_mb=config.get("vcf_addRef", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("vcf_addRef", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("vcf_addRef", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("vcf_addRef", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("vcf_addRef", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("vcf_addRef", {}).get("container", config["default_container"])
    conda:
        "../envs/vcf_addRef.yaml"
    message:
        "{rule}: Add reference to the header of the deepvariant vcf: {input.vcf}"
    script:
        "../scripts/ref_vcf.py"


# rule vcf_changeM2MT:
#     input:
#         "vcf_final/{sample}_{type}_ref.vcf",
#     output:
#         vcf=temp("vcf_final/{sample}_{type}.vcf"),
#     log:
#         "vcf_final/{sample}_{type}_chrMT.log",
#     resources:
#         mem_mb=config.get("vcf_changeM2MT", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
#         mem_per_cpu=config.get("vcf_changeM2MT", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
#         partition=config.get("vcf_changeM2MT", {}).get("partition", config["default_resources"]["partition"]),
#         threads=config.get("vcf_changeM2MT", {}).get("threads", config["default_resources"]["threads"]),
#         time=config.get("vcf_changeM2MT", {}).get("time", config["default_resources"]["time"]),
#     shell:
#         """( awk '{{gsub(/chrM/,"chrMT"); print}}' {input} > {output.vcf} ) &> {log}"""

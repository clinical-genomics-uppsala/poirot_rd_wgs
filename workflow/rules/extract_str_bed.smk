___author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2023"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule extract_str_bed:
    input:
        vcf="cnv_sv/stranger/{sample}_{type}.stranger.vcf",
    output:
        bed=temp("cnv_sv/stranger/{sample}_{type}.stranger.bed"),
    log:
        "cnv_sv/stranger/{sample}_{type}.stranger.bed.log",
    benchmark:
        repeat(
            "cnv_sv/stranger/{sample}_{type}.stranger.bed.benchmark.tsv",
            config.get("extract_str_bed", {}).get("benchmark_repeats", 1),
        )
    resources:
        mem_mb=rule_resource("extract_str_bed", "mem_mb"),
        mem_per_cpu=rule_resource("extract_str_bed", "mem_per_cpu"),
        partition=rule_resource("extract_str_bed", "partition"),
        threads=rule_resource("extract_str_bed", "threads"),
        time=rule_resource("extract_str_bed", "time"),
    container:
        rule_container("extract_str_bed")
    message:
        "{rule}: Convert stranger annotated {input.vcf} to bed file format"
    script:
        "../scripts/extract_str_bed.py"


use rule extract_str_bed as extract_str_bed_panel with:
    input:
        vcf="cnv_sv/stranger/{sample}_{type}.stranger.vcf",
        panel_list=lambda wildcards: get_str_panel_list(wildcards),
    output:
        bed=temp("cnv_sv/stranger/{sample}_{type}_{panel}.stranger.bed"),
    log:
        "cnv_sv/stranger/{sample}_{type}_{panel}.stranger.bed.log",
    benchmark:
        repeat(
            "cnv_sv/stranger/{sample}_{type}_{panel}.stranger.bed.benchmark.tsv",
            config.get("extract_str_bed_panel", {}).get("benchmark_repeats", 1),
        )
    message:
        "{rule}: Convert stranger annotated {input.vcf} to bed file format for loci in {input.panel_list}"

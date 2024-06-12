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
        mem_mb=config.get("extract_str_bed", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("extract_str_bed", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("extract_str_bed", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("extract_str_bed", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("extract_str_bed", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("extract_str_bed", {}).get("container", config["default_container"])
    message:
        "{rule}: Convert stranger annotated {input.vcf} to bed file format"
    script:
        "../scripts/extract_str_bed.py"


use rule extract_str_bed as extract_str_bed_panel with:
    input:
        vcf="cnv_sv/stranger/{sample}_{type}.stranger.vcf",
        panel_list=expand("{panel_dir}/{{panel}}.list", panel_dir=config.get("reference", {}).get("str_panels_dir", "")),
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

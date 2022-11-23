___author__ = "Arielle Reivant Munters, Jessika Nordin"
__copyright__ = "Copyright 2022"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"

rule mosdepth_bedtools:
    input:
        perBase="qc/mosdepth/{sample}_{type}.per-base.bed.gz",
        index="qc/mosdepth/{sample}_{type}.per-base.bed.gz.csi",
        bed=config["reference"]["coverage_bed"],
    output:
        out.=temp("qc/mosdepth_bedtools/{sample}_{type}.mosdepth.lowCov.regions.txt"),
    log:
        "logs/qc/mosdepth_bedtools_{sample}_{type}.log",
    benchmark:
        repeat(
            "qc/mosdepth_bedtools/mosdepth_bedtools_{sample}_{type}.benchmark.tsv",
            config.get("mosdepth_bedtools", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("mosdepth_bedtools", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("mosdepth_bedtools", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("mosdepth_bedtools", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("mosdepth_bedtools", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("mosdepth_bedtools", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("mosdepth_bedtools", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("mosdepth_bedtools", {}).get("container", config["default_container"])
    message:
        "{rule}: Run bedtools to intersect mosdepth with bed of the exons of genes"
    shell:
        "(bedtools intersect -a {input.perBase} -b {input.bed} > {output.out}) &> {log}"


rule create_cov_excel:
    input:
        summary="qc/mosdepth/{sample}_{type}.mosdepth.summary.txt",
        config=config["create_cov_excel"]["covLimits"],
    output:
        out=temp("qc/create_cov_excel/{sample}_{type}.mosdepth.regions.txt"),
    log:
        "logs/qc/create_cov_excel_{sample}_{type}.log",
    benchmark:
        repeat(
            "qc/create_cov_excel/create_cov_excel_{sample}_{type}.benchmark.tsv",
            config.get("create_cov_excel", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("create_cov_excel", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("create_cov_excel", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("create_cov_excel", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("create_cov_excel", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("create_cov_excel", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("create_cov_excel", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("create_cov_excel", {}).get("container", config["default_container"])
    message:
        "{rule}: Get coverage analysis per gene into excel, with tab for each panel and one for all genes in bed"
    script:
        "../scripts/create_excel.py {input.summary} {output.out} {input.config}"

___author__ = "Arielle Reivant Munters, Jessika Nordin"
__copyright__ = "Copyright 2022"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"


rule mosdepth_bedtools:
    input:
        perBase="qc/mosdepth_bed/{sample}_{type}.per-base.bed.gz",
        bed=config["reference"]["coverage_bed"],
    output:
        out=temp("qc/mosdepth_bed/{sample}_{type}.mosdepth.lowCov.regions.txt"),
    log:
        "qc/mosdepth_bed/{sample}_{type}.lowCov.log",
    benchmark:
        repeat(
            "qc/mosdepth_bedtools/mosdepth_bedtools_{sample}_{type}.benchmark.tsv",
            config.get("mosdepth_bedtools", {}).get("benchmark_repeats", 1),
        )
    threads: rule_resource("mosdepth_bedtools", "threads")
    resources:
        mem_mb=rule_resource("mosdepth_bedtools", "mem_mb"),
        mem_per_cpu=rule_resource("mosdepth_bedtools", "mem_per_cpu"),
        partition=rule_resource("mosdepth_bedtools", "partition"),
        threads=rule_resource("mosdepth_bedtools", "threads"),
        time=rule_resource("mosdepth_bedtools", "time"),
    container:
        rule_container("mosdepth_bedtools")
    message:
        "{rule}: Run bedtools to intersect mosdepth with bed of the exons of genes"
    shell:
        "(bedtools intersect -a {input.perBase} -b {input.bed} > {output.out}) &> {log}"


rule create_cov_excel:
    input:
        bedfile=config["reference"]["coverage_bed"],
        cov_regions="qc/mosdepth_bed/{sample}_{type}.regions.bed.gz",
        cov_thresh="qc/mosdepth_bed/{sample}_{type}.thresholds.bed.gz",
        duplication_file="qc/picard_collect_duplication_metrics/{sample}_{type}.duplication_metrics.txt",
        genepanels=config["reference"]["genepanels"],
        low_cov="qc/mosdepth_bed/{sample}_{type}.mosdepth.lowCov.regions.txt",
        summary="qc/mosdepth_bed/{sample}_{type}.mosdepth.summary.txt",
    output:
        out=temp("qc/create_cov_excel/{sample}_{type}.coverage.xlsx"),
    log:
        "qc/create_cov_excel/{sample}_{type}.log",
    benchmark:
        repeat(
            "qc/create_cov_excel/create_cov_excel_{sample}_{type}.benchmark.tsv",
            config.get("create_cov_excel", {}).get("benchmark_repeats", 1),
        )
    threads: rule_resource("create_cov_excel", "threads")
    resources:
        mem_mb=rule_resource("create_cov_excel", "mem_mb"),
        mem_per_cpu=rule_resource("create_cov_excel", "mem_per_cpu"),
        partition=rule_resource("create_cov_excel", "partition"),
        threads=rule_resource("create_cov_excel", "threads"),
        time=rule_resource("create_cov_excel", "time"),
    container:
        rule_container("create_cov_excel")
    message:
        "{rule}: Get coverage analysis per gene into excel, with tab for each panel and one for all genes in bed"
    script:
        "../scripts/create_excel.py"

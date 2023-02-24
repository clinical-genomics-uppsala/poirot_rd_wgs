__author__ = "Arielle R Munters, Jessika Nordin"
__copyright__ = "Copyright 2022, Arielle R Munters"
__email__ = "arielle.munters@scilifelab.u.se"
__license__ = "GPL-3"


rule cp_vcf_all:
    input:
        "vcf_final/{sample}_N.vcf.gz",
    output:
        "results/{sample}/{sample}_snv_indels.vcf.gz",
    params:
        extra=config.get("cp_results", {}).get("extra", ""),
    log:
        "results/logs/{sample}.snv_indels.log",
    benchmark:
        repeat(
            "results/benchmark/{sample}.snv_indels.benchmark.tsv",
            config.get("cp_vcf_all", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cp_results", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_results", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_results", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_results", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_results", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_results", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_results", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    message:
        "{rule}: Move snv_indel final vcf file to result folder for transfer to hospital"
    shell:
        "cp {input} {output}"


rule cp_tbi_all:
    input:
        "vcf_final/{sample}_N.vcf.gz.tbi",
    output:
        "results/{sample}/{sample}_snv_indels.vcf.gz.tbi",
    params:
        extra=config.get("cp_tbi_all", {}).get("extra", ""),
    log:
        "results/logs/{sample}.snv_indels_index.log",
    benchmark:
        repeat(
            "results/benchmark/{sample}.snv_indels_index.benchmark.tsv",
            config.get("cp_results", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cp_results", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_results", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_results", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_results", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_results", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_results", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_results", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    message:
        "{rule}: Move index for snv_indel final vcf file to result folder for transfer to hospital"
    shell:
        "cp {input} {output}"


rule cp_cram_all:
    input:
        "compression/crumble/{sample}_N.crumble.cram",
    output:
        "results/{sample}/{sample}.crumble.cram",
    params:
        extra=config.get("cp_cram_all", {}).get("extra", ""),
    log:
        "results/logs/{sample}.cram.log",
    benchmark:
        repeat(
            "results/benchmark/{sample}.cram.benchmark.tsv",
            config.get("cp_results", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cp_results", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_results", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_results", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_results", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_results", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_results", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_results", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    message:
        "{rule}: Move crumble cram file to result folder for transfer to hospital"
    shell:
        "cp {input} {output}"


rule cp_crai_all:
    input:
        "compression/crumble/{sample}_N.crumble.cram.crai",
    output:
        "results/{sample}/{sample}.crumble.cram.crai",
    params:
        extra=config.get("cp_cram_all", {}).get("extra", ""),
    log:
        "results/logs/{sample}.cram_index.log",
    benchmark:
        repeat(
            "results/benchmark/{sample}.cram_index.benchmark.tsv",
            config.get("cp_results", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cp_results", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_results", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_results", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_results", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_results", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_results", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_results", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    message:
        "{rule}: Move crumble cram index file to result folder for transfer to hospital"
    shell:
        "cp {input} {output}"


rule cp_multiqc:
    input:
        "qc/multiqc/multiqc_DNA.html",
    output:
        "results/multiqc_DNA.html",
    params:
        extra=config.get("cp_multiqc", {}).get("extra", ""),
    log:
        "results/logs/multiqc.log",
    benchmark:
        repeat(
            "results/benchmark/multiqc.benchmark.tsv",
            config.get("cp_results", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cp_results", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_results", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_results", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_results", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_results", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_results", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_results", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    message:
        "{rule}: Move multiqc report to result folder for transfer to hospital"
    shell:
        "cp {input} {output}"


rule cp_spring:
    input:
        "compression/spring/{sample}_{flowcell}_{lane}_{barcode}_{type}.spring",
    output:
        "results/{sample}/spring/{sample}_{flowcell}_{lane}_{barcode}_{type}.spring",
    params:
        extra=config.get("cp_multiqc", {}).get("extra", ""),
    log:
        "results/logs/{sample}_{flowcell}_{lane}_{barcode}_{type}.spring.log",
    benchmark:
        repeat(
            "results/benchmark/{sample}_{flowcell}_{lane}_{barcode}_{type}.spring.benchmark.tsv",
            config.get("cp_results", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cp_results", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_results", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_results", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_results", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_results", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_results", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_results", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    message:
        "{rule}: Move spring files to a spring folder in the result folder for transfer to hospital"
    shell:
        "cp {input} {output}"


rule cp_smn_json:
    input:
        "cnv_sv/smn_caller/{sample}_N.json",
    output:
        "results/{sample}/SMNCopyNumberCaller/{sample}.smn_caller.json",
    params:
        extra=config.get("cp_smn_json", {}).get("extra", ""),
    log:
        "results/logs/{sample}.smn_json.log",
    benchmark:
        repeat(
            "results/benchmark/{sample}.smn_json.benchmark.tsv",
            config.get("cp_results", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cp_results", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_results", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_results", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_results", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_results", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_results", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_results", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    message:
        "{rule}: Move smn_json files to a SMNCopyNumberCaller folder in each sample folder in the result folder for transfer to hospital"
    shell:
        "cp {input} {output}"


rule cp_smn_tsv:
    input:
        "cnv_sv/smn_caller/{sample}_N.tsv",
    output:
        "results/{sample}/SMNCopyNumberCaller/{sample}.smn_caller.tsv",
    params:
        extra=config.get("cp_smn_tsv", {}).get("extra", ""),
    log:
        "results/logs/{sample}.smn_tsv.log",
    benchmark:
        repeat(
            "results/benchmark/{sample}.smn_tsv.benchmark.tsv",
            config.get("cp_results", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cp_results", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_results", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_results", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_results", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_results", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_results", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_results", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    message:
        "{rule}:  Move smn_tsv files to a SMNCopyNumberCaller folder in each sample folder in the result folder for transfer to hospital"
    shell:
        "cp {input} {output}"


rule cp_smn_pdf:
    input:
        "cnv_sv/smn_charts/smn_{sample}_N.pdf",
    output:
        "results/{sample}/SMNCopyNumberCaller/{sample}.smn_charts.pdf",
    params:
        extra=config.get("cp_smn_pdf", {}).get("extra", ""),
    log:
        "results/logs/{sample}.smn_pdf.log",
    benchmark:
        repeat(
            "results/benchmark/{sample}.smn_pdf.benchmark.tsv",
            config.get("cp_results", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cp_results", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_results", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_results", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_results", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_results", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_results", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_results", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    message:
        "{rule}: Move smn_pdf files to a SMNCopyNumberCaller folder in each sample folder in the result folder for transfer to hospital"
    shell:
        "cp {input} {output}"


rule cp_coverage:
    input:
        "qc/create_cov_excel/{sample}_N.coverage.xlsx",
    output:
        "results/{sample}/{sample}.coverage_analysis.xlsx",
    params:
        extra=config.get("cp_coverage", {}).get("extra", ""),
    log:
        "results/logs/{sample}.coverage.log",
    benchmark:
        repeat(
            "results/benchmark/{sample}.coverage.benchmark.tsv",
            config.get("cp_results", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cp_results", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_results", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_results", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_results", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_results", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_results", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_results", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    message:
        "{rule}: Move coverage analysis reports to a sample folder in the result folder for transfer to hospital"
    shell:
        "cp {input} {output}"


rule cp_stranger:
    input:
        "cnv_sv/stranger/{sample}_N.stranger.vcf.gz",
    output:
        "results/{sample}/{sample}.expansionhunter_stranger.vcf.gz",
    params:
        extra=config.get("cp_stranger", {}).get("extra", ""),
    log:
        "results/logs/{sample}.stranger.log",
    benchmark:
        repeat(
            "results/benchmark/{sample}.stranger.benchmark.tsv",
            config.get("cp_results", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cp_results", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_results", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_results", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_results", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_results", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_results", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_results", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    message:
        "{rule}: Move stranger files to a sample folder in the result folder for transfer to hospital"
    shell:
        "cp {input} {output}"


rule cp_reviewer:
    input:
        "cnv_sv/reviewer/{sample}_N/",
    output:
        directory("results/{sample}/expansionhunter_reviewer/"),
    params:
        extra=config.get("cp_reviewer", {}).get("extra", ""),
    log:
        "results/logs/{sample}.reviewer.log",
    benchmark:
        repeat(
            "results/benchmark/{sample}.reviewer.benchmark.tsv",
            config.get("cp_results", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cp_results", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_results", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_results", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_results", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_results", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_results", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_results", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    message:
        "{rule}: Move reviewer folder to a sample folder in the result folder for transfer to hospital"
    shell:
        "cp -r {input} {output}"


rule cp_contamination:
    input:
        "mitochondrial/haplocheck/{sample}_N.contamination.html",
    output:
        "results/{sample}/{sample}.contamination.html",
    params:
        extra=config.get("cp_contamination", {}).get("extra", ""),
    log:
        "results/logs/{sample}.contamination.log",
    benchmark:
        repeat(
            "results/benchmark/{sample}.contamination.benchmark.tsv",
            config.get("cp_results", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cp_results", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_results", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_results", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_results", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_results", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_results", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_results", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    message:
        "{rule}: Move contamination files to a sample folder in the result folder for transfer to hospital"
    shell:
        "cp {input} {output}"


rule cp_merge:
    input:
        "cnv_sv/svdb_merge/{sample}_N.merged.vcf.gz",
    output:
        "results/{sample}/cnv_sv/{sample}.svdb_merged.vcf.gz",
    params:
        extra=config.get("cp_merge", {}).get("extra", ""),
    log:
        "results/logs/{sample}.merge.log",
    benchmark:
        repeat(
            "results/benchmark/{sample}.merge.benchmark.tsv",
            config.get("cp_results", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cp_results", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_results", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_results", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_results", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_results", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_results", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_results", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    message:
        "{rule}: Move SVDB merged files from manta, tiddit and cnvpytor to a cnv_sv folder for each sample in the result folder for transfer to hospital"
    shell:
        "cp {input} {output}"


rule cp_manta:
    input:
        "cnv_sv/manta_run_workflow_n/{sample}/results/variants/diploidSV.vcf.gz",
    output:
        "results/{sample}/cnv_sv/{sample}.manta_diploidSV.vcf.gz",
    params:
        extra=config.get("cp_manta", {}).get("extra", ""),
    log:
        "results/logs/{sample}.manta.log",
    benchmark:
        repeat(
            "results/benchmark/{sample}.manta.benchmark.tsv",
            config.get("cp_results", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cp_results", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_results", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_results", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_results", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_results", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_results", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_results", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    message:
        "{rule}: Move manta files to a cnv_sv folder for each sample in the result folder for transfer to hospital"
    shell:
        "cp {input} {output}"


rule cp_tiddit:
    input:
        "cnv_sv/tiddit/{sample}_N.vcf.gz",
    output:
        "results/{sample}/cnv_sv/{sample}.tiddit.vcf.gz",
    params:
        extra=config.get("cp_tiddit", {}).get("extra", ""),
    log:
        "results/logs/{sample}.tiddit.log",
    benchmark:
        repeat(
            "results/benchmark/{sample}.tiddit.benchmark.tsv",
            config.get("cp_results", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cp_results", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_results", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_results", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_results", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_results", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_results", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_results", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    message:
        "{rule}: Move tiddit files to a cnv_sv folder for each sample in the result folder for transfer to hospital"
    shell:
        "cp {input} {output}"


rule cp_cnvpytor:
    input:
        "cnv_sv/cnvpytor/{sample}_N.vcf.gz",
    output:
        "results/{sample}/cnv_sv/{sample}.cnvpytor.vcf.gz",
    params:
        extra=config.get("cp_cnvpytor", {}).get("extra", ""),
    log:
        "results/logs/{sample}.cnvpytor.log",
    benchmark:
        repeat(
            "results/benchmark/{sample}.cnvpytor.benchmark.tsv",
            config.get("cp_results", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cp_results", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_results", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_results", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_results", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_results", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_results", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_results", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    message:
        "{rule}: Move cnvpytor files to a cnv_sv folder for each sample in the result folder for transfer to hospital"
    shell:
        "cp {input} {output}"


rule cp_cnvpytor_filter:
    input:
        "cnv_sv/cnvpytor/{sample}_N.filtered.vcf.gz",
    output:
        "results/{sample}/cnv_sv/{sample}.cnvpytor_filtered.vcf.gz",
    params:
        extra=config.get("cp_cnvpytor_filter", {}).get("extra", ""),
    log:
        "results/logs/{sample}.cnvpytor_filter.log",
    benchmark:
        repeat(
            "results/benchmark/{sample}.cnvpytor_filter.benchmark.tsv",
            config.get("cp_results", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("cp_results", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_results", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_results", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_results", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_results", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_results", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_results", {}).get("container", config["default_container"])
    conda:
        "../envs/cp_results.yaml"
    message:
        "{rule}: Move filtered cnvpytor files to a cnv_sv folder for each sample in the result folder for transfer to hospital"
    shell:
        "cp {input} {output}"

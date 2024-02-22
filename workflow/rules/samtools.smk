


rule samtools_extract_non_chr_reads:
    input:
        bam="alignment/bwa_mem/{sample}_{type}.bam",
        bai="alignment/bwa_mem/{sample}_{type}.bam.bai",
    output:
        bam=temp("alignment/samtools_extract_reads/{sample}_{type}_non_chr.bam"),
    params:
        contigs=non_chr_contigs,
        extra=config.get("samtools_extract_reads", {}).get("extra", ""),
    log:
        "alignment/samtools_extract_reads/{sample}_{type}_non_chr.bam.log",
    benchmark:
        repeat(
            "alignment/samtools_extract_reads/{sample}_{type}_non_chr.bam.benchmark.tsv",
            config.get("samtools_extract_reads", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("samtools_extract_reads", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("samtools_extract_reads", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("samtools_extract_reads", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("samtools_extract_reads", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("samtools_extract_reads", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("samtools_extract_reads", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("samtools_extract_reads", {}).get("container", config["default_container"])
    message:
        "{rule}: create bam {output} with only reads from {params.contigs}"
    shell:
        "(samtools view -@ {threads} {params.extra} -b {input} {params.contigs} > {output}) &> {log}"
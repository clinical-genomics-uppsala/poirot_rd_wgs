__author__ = "Arielle R Munters, Jessika Nordin"
__copyright__ = "Copyright 2022, Arielle R Munters"
__email__ = "arielle.munters@scilifelab.u.se"
__license__ = "GPL-3"


rule cp_vcf_all:
    input:
        "vcf_final/{sample}.vcf.gz",
    output:
        "results/{sample}/{sample}_snv_indels.vcf.gz",
    params:
        extra=config.get("cp_vcf_all", {}).get("extra", ""),
    log:
        "results/{sample}.snv_indels.log",
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



mv vcf_final/{sample}.vcf.gz results/{sample}/

samtools index compression/crumble/{sample}_{type}.crumble.cram

mv compression/crumble/{sample}_{type}.crumble.cram results/{sample}/
mv compression/spring/{sample}_{flowcell}_{lane}_{barcode}_{type}.spring results/{sample}/
mv cnv_sv/svdb_merge/{sample}_{type}.merged.vcf results/{sample}/{sample}_{type}.svdb_merged.vcf
mv cnv_sv/stranger/{sample}_{type}.stranger.vcf results/{sample}/{sample}_{type}.expansionhunter_stranger.vcf
mv cnv_sv/reviewer/{sample} results/{sample}/expansionhunter_reviewer
mv mitochondrial/haplocheck/{sample}_{type}.contamination.html results/{sample}/
mv qc/create_cov_excel/{sample}_{type}.coverage.xlsx results/{sample}/{sample}_{type}.coverage_analysis.xlsx
mv cnv_sv/smn_caller/{sample}_{type}.json results/{sample}/{sample}_{type}.smn_caller.json
mv cnv_sv/smn_caller/{sample}_{type}.tsv results/{sample}/{sample}_{type}.smn_caller.tsv
mv cnv_sv/smn_charts/{sample}_{type}.pdf results/{sample}/{sample}_{type}.smn_charts.pdf
mv qc/multiqc/multiqc_DNA.html results/





rule cp_tbi_all:
    input:
        "vcf_final/{sample}.vcf.gz.tbi",
    output:
        "results/{sample}/{sample}_snv_indels.vcf.gz.tbi",
    params:
        extra=config.get("cp_tbi_all", {}).get("extra", ""),
    log:
        "results/{sample}.snv_indels_index.log",
    benchmark:
        repeat(
            "results/benchmark/{sample}.snv_indels_index.benchmark.tsv",
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
        "{rule}: Move index for snv_indel final vcf file to result folder for transfer to hospital"
    shell:
        "cp {input} {output}"


rule cp_vcf_aml:
    input:
        "parabricks/pbrun_mutectcaller_tn/{sample}.vep.aml.vcf.gz",
    output:
        "Results/{project}/{sample}/SNV_indels/{sample}_TN.vep.aml.vcf.gz",
    threads: config.get("cp_vcf_aml", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_vcf_aml", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_vcf_aml", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_vcf_aml", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_vcf_aml", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_vcf_aml", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_vcf_aml", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_tbi_aml:
    input:
        "parabricks/pbrun_mutectcaller_tn/{sample}.vep.aml.vcf.gz.tbi",
    output:
        "Results/{project}/{sample}/SNV_indels/{sample}_TN.vep.aml.vcf.gz.tbi",
    threads: config.get("cp_tbi_aml", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_tbi_aml", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_tbi_aml", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_tbi_aml", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_tbi_aml", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_tbi_aml", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_tbi_aml", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_vcf_tn:
    input:
        "parabricks/pbrun_mutectcaller_tn/{sample}.vep.vcf.gz",
    output:
        "Results/{project}/{sample}/SNV_indels/{sample}_TN.vep.vcf.gz",
    threads: config.get("cp_vcf_tn", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_vcf_tn", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_vcf_tn", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_vcf_tn", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_vcf_tn", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_vcf_tn", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_vcf_tn", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_tbi_tn:
    input:
        "parabricks/pbrun_mutectcaller_tn/{sample}.vep.vcf.gz.tbi",
    output:
        "Results/{project}/{sample}/SNV_indels/{sample}_TN.vep.vcf.gz.tbi",
    threads: config.get("cp_tbi_tn", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_tbi_tn", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_tbi_tn", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_tbi_tn", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_tbi_tn", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_tbi_tn", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_tbi_tn", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_vcf_t:
    input:
        "parabricks/pbrun_mutectcaller_t/{sample}_T.vep.vcf.gz",
    output:
        "Results/{project}/{sample}/SNV_indels/{sample}_T.vep.vcf.gz",
    threads: config.get("cp_vcf_t", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_vcf_t", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_vcf_t", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_vcf_t", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_vcf_t", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_vcf_t", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_vcf_t", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_tbi_t:
    input:
        "parabricks/pbrun_mutectcaller_t/{sample}_T.vep.vcf.gz.tbi",
    output:
        "Results/{project}/{sample}/SNV_indels/{sample}_T.vep.vcf.gz.tbi",
    threads: config.get("cp_tbi_t", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_tbi_t", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_tbi_t", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_tbi_t", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_tbi_t", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_tbi_t", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_tbi_t", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_cram:
    input:
        "compression/crumble/{sample}_{type}.crumble.cram",
    output:
        "Results/{project}/{sample}/Cram/{sample}_{type}.crumble.cram",
    threads: config.get("cp_cram", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_cram", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_cram", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_cram", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_cram", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_cram", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_cram", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_crai:
    input:
        "compression/crumble/{sample}_{type}.crumble.cram.crai",
    output:
        "Results/{project}/{sample}/Cram/{sample}_{type}.crumble.cram.crai",
    threads: config.get("cp_crai", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_crai", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_crai", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_crai", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_crai", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_crai", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_crai", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_tsv_mutectcaller_all:
    input:
        "tsv_files/{sample}_mutectcaller_tn.all.tsv",
    output:
        "Results/{project}/{sample}/SNV_indels/{sample}_mutectcaller_TN.all.tsv",
    threads: config.get("cp_tsv_mutectcaller_all", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_tsv_mutectcaller_all", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_tsv_mutectcaller_all", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_tsv_mutectcaller_all", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_tsv_mutectcaller_all", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_tsv_mutectcaller_all", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_tsv_mutectcaller_all", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_tsv_manta_all:
    input:
        "tsv_files/{sample}_manta_tn.all.tsv",
    output:
        "Results/{project}/{sample}/SV/{sample}_manta_TN.all.tsv",
    threads: config.get("cp_tsv_manta_all", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_tsv_manta_all", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_tsv_manta_all", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_tsv_manta_all", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_tsv_manta_all", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_tsv_manta_all", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_tsv_manta_all", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_tsv_mutectcaller_aml:
    input:
        "tsv_files/{sample}_mutectcaller_tn.aml.tsv",
    output:
        "Results/{project}/{sample}/SNV_indels/{sample}_mutectcaller_TN.aml.tsv",
    threads: config.get("cp_tsv_muntectcaller_aml", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_tsv_muntectcaller_aml", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_tsv_muntectcaller_aml", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_tsv_muntectcaller_aml", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_tsv_muntectcaller_aml", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_tsv_muntectcaller_aml", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_tsv_muntectcaller_aml", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_tsv_manta_aml:
    input:
        "tsv_files/{sample}_manta_tn.aml.tsv",
    output:
        "Results/{project}/{sample}/SV/{sample}_manta_TN.aml.tsv",
    threads: config.get("cp_tsv_manta_aml", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_tsv_manta_aml", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_tsv_manta_aml", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_tsv_manta_aml", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_tsv_manta_aml", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_tsv_manta_aml", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_tsv_manta_aml", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_manta_vcf:
    input:
        "cnv_sv/manta_run_workflow_tn/{sample}.ssa.vcf.gz",
    output:
        "Results/{project}/{sample}/SV/{sample}_TN.ssa.vcf.gz",
    threads: config.get("cp_manta_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_manta_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_manta_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_manta_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_manta_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_manta_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_manta_vcf", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_manta_tbi:
    input:
        "cnv_sv/manta_run_workflow_tn/{sample}.ssa.vcf.gz.tbi",
    output:
        "Results/{project}/{sample}/SV/{sample}_TN.ssa.vcf.gz.tbi",
    threads: config.get("cp_manta_tbi", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_manta_tbi", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_manta_tbi", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_manta_tbi", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_manta_tbi", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_manta_tbi", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_manta_tbi", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_manta_all_vcf:
    input:
        "cnv_sv/manta_run_workflow_tn/{sample}.ssa.all.vcf.gz",
    output:
        "Results/{project}/{sample}/SV/{sample}_TN.ssa.all.vcf.gz",
    threads: config.get("cp_manta_all_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_manta_all_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_manta_all_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_manta_all_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_manta_all_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_manta_all_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_manta_all_vcf", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_manta_all_tbi:
    input:
        "cnv_sv/manta_run_workflow_tn/{sample}.ssa.all.vcf.gz.tbi",
    output:
        "Results/{project}/{sample}/SV/{sample}_TN.ssa.all.vcf.gz.tbi",
    threads: config.get("cp_manta_all_tbi", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_manta_all_tbi", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_manta_all_tbi", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_manta_all_tbi", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_manta_all_tbi", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_manta_all_tbi", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_manta_all_tbi", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_manta_aml_vcf:
    input:
        "cnv_sv/manta_run_workflow_tn/{sample}.ssa.aml.vcf.gz",
    output:
        "Results/{project}/{sample}/SV/{sample}_TN.ssa.aml.vcf.gz",
    threads: config.get("cp_manta_aml_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_manta_aml_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_manta_aml_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_manta_aml_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_manta_aml_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_manta_aml_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_manta_aml_vcf", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_manta_aml_tbi:
    input:
        "cnv_sv/manta_run_workflow_tn/{sample}.ssa.aml.vcf.gz.tbi",
    output:
        "Results/{project}/{sample}/SV/{sample}_TN.ssa.aml.vcf.gz.tbi",
    threads: config.get("cp_manta_aml_tbi", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_manta_aml_tbi", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_manta_aml_tbi", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_manta_aml_tbi", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_manta_aml_tbi", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_manta_aml_tbi", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_manta_aml_tbi", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_pindel_vcf:
    input:
        "cnv_sv/pindel/{sample}.vcf.gz",
    output:
        "Results/{project}/{sample}/SNV_indels/{sample}.pindel.vcf.gz",
    threads: config.get("cp_pindel_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_pindel_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_pindel_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_pindel_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_pindel_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_pindel_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_pindel_vcf", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_pindel_tbi:
    input:
        "cnv_sv/pindel/{sample}.vcf.gz.tbi",
    output:
        "Results/{project}/{sample}/SNV_indels/{sample}.pindel.vcf.gz.tbi",
    threads: config.get("cp_pindel_tbi", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_pindel_tbi", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_pindel_tbi", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_pindel_tbi", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_pindel_tbi", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_pindel_tbi", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_pindel_tbi", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_cnvkit_table:
    input:
        "cnv_sv/cnvkit_table/{sample}_T.CNV.xlsx",
    output:
        "Results/{project}/{sample}/CNV/{sample}_{type}.CNV.xlsx",
    threads: config.get("cp_cnvkit_table", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_cnvkit_table", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_cnvkit_table", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_cnvkit_table", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_cnvkit_table", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_cnvkit_table", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_cnvkit_table", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_cnvkit_vcf:
    input:
        "cnv_sv/cnvkit_vcf/{sample}_T.vcf.gz",
    output:
        "Results/{project}/{sample}/CNV/{sample}_T.vcf.gz",
    threads: config.get("cp_cnvkit_vcf", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_cnvkit_vcf", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_cnvkit_vcf", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_cnvkit_vcf", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_cnvkit_vcf", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_cnvkit_vcf", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_cnvkit_vcf", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_cnvkit_tbi:
    input:
        "cnv_sv/cnvkit_vcf/{sample}_T.vcf.gz.tbi",
    output:
        "Results/{project}/{sample}/CNV/{sample}_T.vcf.gz.tbi",
    threads: config.get("cp_cnvkit_tbi", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_cnvkit_tbi", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_cnvkit_tbi", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_cnvkit_tbi", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_cnvkit_tbi", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_cnvkit_tbi", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_cnvkit_tbi", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_cnvkit_diagram:
    input:
        "cnv_sv/cnvkit_diagram/{sample}_T.png",
    output:
        "Results/{project}/{sample}/CNV/{sample}_T.png",
    threads: config.get("cp_cnvkit_diagram", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_cnvkit_diagram", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_cnvkit_diagram", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_cnvkit_diagram", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_cnvkit_diagram", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_cnvkit_diagram", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_cnvkit_diagram", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_cnvkit_scatter:
    input:
        "cnv_sv/cnvkit_scatter/{sample}_T_chr{chr}.png",
    output:
        "Results/{project}/{sample}/CNV/{sample}_T_chr{chr}.png",
    threads: config.get("cp_cnvkit_scatter", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_cnvkit_scatter", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_cnvkit_scatter", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_cnvkit_scatter", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_cnvkit_scatter", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_cnvkit_scatter", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_cnvkit_scatter", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_multiqc:
    input:
        "qc/multiqc/multiqc_DNA.html",
    output:
        "Results/MultiQC_TN.html",
    threads: config.get("cp_multiqc", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_multiqc", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_multiqc", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_multiqc", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_multiqc", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_multiqc", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_multiqc", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"


rule cp_spring_archive:
    input:
        "compression/spring/{sample}_{flowcell}_{lane}_{barcode}_{type}.spring",
    output:
        "Archive/{project}/{sample}_{flowcell}_{lane}_{barcode}_{type}.spring",
    threads: config.get("cp_spring_archive", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("cp_spring_archive", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("cp_spring_archive", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("cp_spring_archive", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("cp_spring_archive", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("cp_spring_archive", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("cp_spring_archive", {}).get("container", config["default_container"])
    shell:
        "cp {input} {output}"

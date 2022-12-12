

rule glnexus:
    input:
        # gvcf=expand("parabricks/pbrun_deepvariant/{sample}.g.vcf", sample=get_samples(samples)),
        # gvcf=expand("qc/peddy/{sample}_{type}.g.vcf", sample=get_samples(samples)),
        gvcf=[
            "snv_indels/deepvariant/{}_{}.g.vcf".format(sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ],
    output:
        bcf=temp("qc/peddy/all.bcf"),
        glnexus=temp(directory("GLnexus.DB")),
    params:
        extra=config.get("glnexus", {}).get("extra", ""),
        in_gvcf=get_in_gvcf,
    log:
        "qc/peddy/all.bcf.log",
    benchmark:
        repeat(
            "qc/peddy/all.bcf.benchmark.tsv",
            config.get("glnexus", {}).get("benchmark_repeats", 1),
        )
    threads: config.get("glnexus", {}).get("threads", config["default_resources"]["threads"])
    resources:
        mem_mb=config.get("glnexus", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("glnexus", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("glnexus", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("glnexus", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("glnexus", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("glnexus", {}).get("container", config["default_container"])
    message:
        "{rule}: Run GLNexus for joint genotyping of Deepvariant gVCFs"
    shell:
        "glnexus_cli --config DeepVariant {params.extra} -i {params.in_gvcf} > {output.bcf}"


rule bcftools_view:
    input:
        "qc/peddy/all.bcf",
    output:
        temp("qc/peddy/all.vcf.gz"),
    log:
        "qc/peddy/all.vcf.log",
    benchmark:
        repeat(
            "qc/peddy/all.bcf.benchmark.tsv",
            config.get("bcftools_view", {}).get("benchmark_repeats", 1),
        )
    resources:
        mem_mb=config.get("bcftools_view", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("bcftools_view", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("bcftools_view", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("bcftools_view", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("bcftools_view", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("bcftools_view", {}).get("container", config["default_container"])
    message:
        "{rule}: Run bcftools view to convert glNexus bcf to vcf and then bgzip and tabix"
    shell:
        """
        bcftools view -O z -o {output} {input}
        """


rule create_ped:
    input:
        config["peddy"]["samples"],
    output:
        temp("qc/peddy/all.ped"),
    log:
        "qc/peddy/all.ped.log",
    resources:
        mem_mb=config.get("create_ped", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("create_ped", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("create_ped", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("create_ped", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("create_ped", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("create_ped", {}).get("container", config["default_container"])
    message:
        "{rule}: Create a peddy ped/FAM file from the SampleSheet.csv file"
    script:
        "../scripts/create_peddy_fam.py"


rule create_peddy_mqc_tsv:
    input:
        peddy_rel_check="qc/peddy/peddy.ped_check.csv",
        peddy_sex_check="qc/peddy/peddy.sex_check.csv",
        ped="qc/peddy/all.ped",
    output:
        rel_check_mqc=temp("qc/peddy/peddy_rel_check_mqc.tsv"),
        sex_check_mqc=temp("qc/peddy/peddy_sex_check_mqc.tsv"),
    params:
        pre="qc/peddy/peddy_mqc",
    log:
        "qc/peddy/peddy.log",
    benchmark:
        repeat(
            "qc/peddy/create_peddy_mqc_tsv.benchmark.tsv",
            config.get("create_peddy_mqc_tsv", {}).get("benchmark_repeats", 1),
        )
    resources:
        mem_mb=config.get("create_peddy_mqc_tsv", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("create_peddy_mqc_tsv", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("create_peddy_mqc_tsv", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("create_peddy_mqc_tsv", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("create_peddy_mqc_tsv", {}).get("time", config["default_resources"]["time"]),
    message:
        "{rule}: Create multiqc custom content embedded config tsv files from peddy sex_check and ped_check files"
    script:
        "../scripts/create_peddy_mqc_config.py"

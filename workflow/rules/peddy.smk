___author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2022"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule create_ped:
    input:
        config["samples"],
    output:
        temp("qc/peddy/all.ped"),
    log:
        "qc/peddy/all.ped.log",
    resources:
        mem_mb=rule_resource("create_ped", "mem_mb"),
        mem_per_cpu=rule_resource("create_ped", "mem_per_cpu"),
        partition=rule_resource("create_ped", "partition"),
        threads=rule_resource("create_ped", "threads"),
        time=rule_resource("create_ped", "time"),
    container:
        rule_container("create_ped")
    message:
        "{rule}: Create a peddy ped/FAM file from the samples.tsv file"
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
        mem_mb=rule_resource("create_peddy_mqc_tsv", "mem_mb"),
        mem_per_cpu=rule_resource("create_peddy_mqc_tsv", "mem_per_cpu"),
        partition=rule_resource("create_peddy_mqc_tsv", "partition"),
        threads=rule_resource("create_peddy_mqc_tsv", "threads"),
        time=rule_resource("create_peddy_mqc_tsv", "time"),
    container:
        rule_container("create_peddy_mqc_tsv")
    message:
        "{rule}: Create multiqc custom content embedded config tsv files from peddy sex_check and ped_check files"
    script:
        "../scripts/create_peddy_mqc_config.py"

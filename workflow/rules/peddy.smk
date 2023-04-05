___author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2022"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


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
    container:
        config.get("create_peddy_mqc_tsv", {}).get("container", config["default_container"])
    message:
        "{rule}: Create multiqc custom content embedded config tsv files from peddy sex_check and ped_check files"
    script:
        "../scripts/create_peddy_mqc_config.py"

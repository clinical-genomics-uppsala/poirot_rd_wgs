___author__ = "Padraic Corcoran"
__copyright__ = "Copyright 2025"
__email__ = "padraic.corcoran@scilifelab.uu.se"
__license__ = "GPL-3"


rule create_somalier_mqc_tsv:
    input:
        pairs="qc/somalier_trio/somalier_relate.pairs.tsv",
        samples="qc/somalier_trio/somalier_relate.samples.tsv",
        ped="qc/somalier_trio/somalier_all.ped",
    output:
        rel_check_mqc=temp("qc/somalier_trio/somalier_rel_check_mqc.tsv"),
        sex_check_mqc=temp("qc/somalier_trio/somalier_sex_check_mqc.tsv"),
        general_stats_mqc=temp("qc/somalier_trio/somalier_general_stats_mqc.tsv"),
    params:
        script=f"{workflow.basedir}/scripts/create_somalier_mqc_config.py",
        mqc_config=config.get("somalier_trio_mqc", {}).get("mqc_config", ""),
    log:
        "qc/somalier_trio/somalier_mqc.log",
    benchmark:
        repeat(
            "qc/somalier_trio/create_somalier_mqc_tsv.benchmark.tsv",
            config.get("create_somalier_mqc_tsv", {}).get("benchmark_repeats", 1),
        )
    resources:
        mem_mb=config.get("create_somalier_mqc_tsv", {}).get("mem_mb", config["default_resources"]["mem_mb"]),
        mem_per_cpu=config.get("create_somalier_mqc_tsv", {}).get("mem_per_cpu", config["default_resources"]["mem_per_cpu"]),
        partition=config.get("create_somalier_mqc_tsv", {}).get("partition", config["default_resources"]["partition"]),
        threads=config.get("create_somalier_mqc_tsv", {}).get("threads", config["default_resources"]["threads"]),
        time=config.get("create_somalier_mqc_tsv", {}).get("time", config["default_resources"]["time"]),
    container:
        config.get("create_somalier_mqc_tsv", {}).get("container", config["default_container"])
    message:
        "{rule}: Create multiqc custom content embedded config tsv files from somalier sex_check and relatedness files"
    shell:
        """
        exec &> {log}
        set -ex
        
        echo "Starting Somalier MultiQC TSV creation"
        echo "Script path: {params.script}"
        ls -l {params.script}
        
        python3 --version
        
        python3 {params.script} \
            --pairs {input.pairs} \
            --samples {input.samples} \
            --ped {input.ped} \
            --config {params.mqc_config} \
            --rel-check-mqc {output.rel_check_mqc} \
            --sex-check-mqc {output.sex_check_mqc} \
            --general-stats-mqc {output.general_stats_mqc}
        """

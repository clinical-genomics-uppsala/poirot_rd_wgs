___author__ = "Jessika Nordin, Padraic Corcoran"
__copyright__ = "Copyright 2022"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"

import pandas
import yaml

from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *
from snakemake.utils import min_version
from snakemake.utils import validate

min_version("6.10")


### Set and validate config file
configfile: "config/config.yaml"


validate(config, schema="../schemas/config.schema.yaml")

config = load_resources(config, config["resources"])
validate(config, schema="../schemas/resources.schema.yaml")


### Read and validate samples file
samples = pandas.read_table(config["samples"], dtype=str).set_index("sample", drop=False)
validate(samples, schema="../schemas/samples.schema.yaml")

if samples[~pandas.isnull(samples.trio_member)].shape[0] % 3 != 0:
        sys.exit('Not all members of the trios are available in the Sample Sheet')

### Read and validate units file
units = (
    pandas.read_table(config["units"], dtype=str)
    .set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")


### Set wildcard constraints
wildcard_constraints:
    barcode="[A-Z]",
    chr="[^_]",
    flowcell="[A-Z0-9]",
    lane="L[0-9]",
    read="fastq[1|2]",
    sample="|".join(get_samples(samples)),
    type="N|T|R",


### Functions
def get_flowcell(units, wildcards):
    flowcells = set([u.flowcell for u in get_units(units, wildcards)])
    if len(flowcells) > 1:
        raise ValueError("Sample type combination from different sequence flowcells")
    return flowcells.pop()


# def get_in_gvcf(wildcards):
#     gvcf_list = [
#         "snv_indels/deepvariant_peddy/{}_{}.g.vcf".format(sample, t)
#         for sample in get_samples(samples)
#         for t in get_unit_types(units, sample)
#     ]
#     return " -i ".join(gvcf_list)


def get_spring_extra(wildcards: snakemake.io.Wildcards):
    extra = config.get("spring", {}).get("extra", "")
    if get_fastq_file(units, wildcards, "fastq1").endswith(".gz"):
        extra = "%s %s" % (extra, "-g")
    return extra


def get_bam_input(wildcards, use_sample_wildcard=True, use_type_wildcard=True, by_chr=False):

    if use_sample_wildcard and use_type_wildcard is True:
        sample_str = "{}_{}".format(wildcards.sample, wildcards.type)
    elif use_sample_wildcard and use_type_wildcard is not True:
        sample_str = "{}_{}".format(wildcards.sample, "N")
    else:
        sample_str = wildcards.file

    aligner = config.get("aligner", None)
    if aligner is None:
        sys.exit("aligner missing from config, valid options: bwa_gpu or bwa_cpu")
    elif aligner == "bwa_gpu":
        bam_input = "parabricks/pbrun_fq2bam/{}.bam".format(sample_str)
    elif aligner == "bwa_cpu":
        if by_chr:  # if a bam for single chromosome is needed
            bam_input = "alignment/picard_mark_duplicates/{}_{}.bam".format(sample_str, wildcards.chr)
        else:
            bam_input = "alignment/samtools_merge_bam/{}.bam".format(sample_str)
    else:
        sys.exit("valid options for aligner are: bwa_gpu or bwa_cpu")

    bai_input = "{}.bai".format(bam_input)

    return (bam_input, bai_input)


def get_vcf_input(wildcards):

    caller = config.get("snp_caller", None)
    if caller is None:
        sys.exit("snp_caller missing from config, valid options: deepvariant_gpu or deepvariant_cpu")
    elif caller == "deepvariant_gpu":
        vcf_input = "parabricks/pbrun_deepvariant/{}_{}.vcf".format(wildcards.sample, wildcards.type)
    elif caller == "deepvariant_cpu":
        vcf_input = "snv_indels/deepvariant/{}_{}.vcf".format(wildcards.sample, wildcards.type)
    else:
        sys.exit("Invalid options for snp_caller, valid options are: deepvariant_gpu or deepvariant_cpu")

    return vcf_input


def get_gvcf_list(wildcards):

    caller = config.get("snp_caller", None)
    if caller is None:
        sys.exit("snp_caller missing from config, valid options: deepvariant_gpu or deepvariant_cpu")
    elif caller == "deepvariant_gpu":
        gvcf_path = "parabricks/pbrun_deepvariant"
    elif caller == "deepvariant_cpu":
        gvcf_path = "snv_indels/deepvariant"
    else:
        sys.exit("Invalid options for snp_caller, valid options are: deepvariant_gpu or deepvariant_cpu")

    gvcf_list = [
        "{}/{}_{}.g.vcf".format(gvcf_path, sample, t)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
    ]
    

    return gvcf_list


def get_in_gvcf(wildcards):
    gvcf_list = get_gvcf_list(wildcards)
    return " -i ".join(gvcf_list)


def get_spring_extra(wildcards: snakemake.io.Wildcards):
    extra = config.get("spring", {}).get("extra", "")
    if get_fastq_file(units, wildcards, "fastq1").endswith(".gz"):
        extra = "%s %s" % (extra, "-g")
    return extra


def get_parent_bams(wildcards):
    aligner = config.get("aligner", None)

    if aligner is None:
        sys.exit("aligner missing from config, valid options: bwa_gpu or bwa_cpu")
    elif aligner == "bwa_gpu":
        bam_path = "parabricks/pbrun_fq2bam"
    elif aligner == "bwa_cpu":
        bam_path = "alignment/samtools_merge_bam"

    proband_sample = samples[samples.index == wildcards.sample]
    trio_id = proband_sample.at[wildcards.sample, "trioid"]

    mother_sample = samples[(samples.trio_member == "mother") & (samples.trioid == trio_id)].index[0]
    mother_bam = "{}/{}_{}.bam".format(bam_path, mother_sample, list(get_unit_types(units, mother_sample))[0])

    father_sample = samples[(samples.trio_member == "father") & (samples.trioid == trio_id)].index[0]
    father_bam = "{}/{}_{}.bam".format(bam_path, father_sample, list(get_unit_types(units, father_sample))[0])

    bam_list = [mother_bam, father_bam]


    return bam_list



def get_glnexus_input(wildcards, input):
   
    gvcf_input =  "-i {}".format(" -i ".join(input.gvcfs))
   
    return gvcf_input


def compile_output_list(wildcards: snakemake.io.Wildcards):
    files = {
        "cnv_sv/cnvpytor": ["vcf.gz"],
        "cnv_sv/expansionhunter": ["vcf.gz"],
        "cnv_sv/smn_caller": ["tsv"],
        "cnv_sv/stranger": ["stranger.vcf.gz"],
        "cnv_sv/svdb_query": ["svdb_query.vcf.gz"],
        "cnv_sv/tiddit": ["vcf.gz"],
        "compression/crumble": ["crumble.cram.crai"],
        "qc/create_cov_excel": ["coverage.xlsx"],
        "mitochondrial/gatk_select_variants_final": ["vcf"],
        "vcf_final": ["vep_annotated.filter.germline.vcf.gz.tbi"],
    }
    output_files = [
        "%s/%s_%s.%s" % (prefix, sample, unit_type, suffix)
        for prefix in files.keys()
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        for suffix in files[prefix]
    ]
    output_files += [
        "cnv_sv/manta_run_workflow_n/%s/results/variants/diploidSV.vcf.gz" % (sample) for sample in get_samples(samples)
    ]
    output_files += [
        "cnv_sv/reviewer/%s_%s/" % (sample, unit_type)
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
    ]
    output_files += [
        "cnv_sv/automap/%s_%s/%s_%s.HomRegions.tsv" % (sample, unit_type, sample, unit_type)
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
    ]
    output_files += [
        "cnv_sv/smn_charts/smn_%s_%s.pdf" % (sample, unit_type)
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
    ]
    output_files += ["qc/multiqc/multiqc_DNA.html"]
    output_files += ["results/multiqc_DNA.html"]
    output_files += [
        "qc/peddy/peddy.peddy.ped",
        "qc/peddy/peddy.ped_check.csv",
        "qc/peddy/peddy.sex_check.csv",
        "qc/peddy/peddy.het_check.csv",
        "qc/peddy/peddy.html",
        "qc/peddy/peddy.vs.html",
        "qc/peddy/peddy.background_pca.json",
    ]

    files = {
        "snv_indels/glnexus": ["vcf.gz"],
        "cnv_sv/upd": ["upd_regions.bed"],
    }
    output_files += [
        "%s/%s_%s.%s" % (prefix, sample, unit_type,suffix)
        for prefix in files.keys()
        for sample in samples[samples.trio_member == 'proband'].index
        for unit_type in get_unit_types(units, sample)
        for suffix in files[prefix]
    ]

   #output_files += ["snv_indels/deeptrio/61743_N/child.g.vcf","snv_indels/deeptrio/61743_N/parent1.g.vcf"]

    output_files += [
        "compression/spring/%s_%s_%s_%s_%s.spring" % (sample, flowcell, lane, barcode, t)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
        for flowcell in set(
            [
                u.flowcell
                for u in units.loc[
                    (
                        sample,
                        t,
                    )
                ]
                .dropna()
                .itertuples()
            ]
        )
        for barcode in set(
            [
                u.barcode
                for u in units.loc[
                    (
                        sample,
                        t,
                    )
                ]
                .dropna()
                .itertuples()
            ]
        )
        for lane in set(
            [
                u.lane
                for u in units.loc[
                    (
                        sample,
                        t,
                    )
                ]
                .dropna()
                .itertuples()
            ]
        )
    ]
    # output_files = [
    #     "results/%s/spring/%s_%s_%s_%s_%s.spring" % (sample, sample, flowcell, lane, barcode, t)
    #     for sample in get_samples(samples)
    #     for t in get_unit_types(units, sample)
    #     for flowcell in set(
    #         [
    #             u.flowcell
    #             for u in units.loc[
    #                 (
    #                     sample,
    #                     t,
    #                 )
    #             ]
    #             .dropna()
    #             .itertuples()
    #         ]
    #     )
    #     for barcode in set(
    #         [
    #             u.barcode
    #             for u in units.loc[
    #                 (
    #                     sample,
    #                     t,
    #                 )
    #             ]
    #             .dropna()
    #             .itertuples()
    #         ]
    #     )
    #     for lane in set(
    #         [
    #             u.lane
    #             for u in units.loc[
    #                 (
    #                     sample,
    #                     t,
    #                 )
    #             ]
    #             .dropna()
    #             .itertuples()
    #         ]
    #     )
    # ]

    # output_files = ["results/%s/cnv_sv/%s.cnvpytor_filtered.vcf.gz" % (sample, sample) for sample in get_samples(samples)]
    # output_files = ["results/%s/cnv_sv/%s.cnvpytor.vcf.gz" % (sample, sample) for sample in get_samples(samples)]
    # output_files = ["results/%s/cnv_sv/%s.tiddit.vcf.gz" % (sample, sample) for sample in get_samples(samples)]
    # output_files = ["results/%s/cnv_sv/%s.manta_diploidSV.vcf.gz" % (sample, sample) for sample in get_samples(samples)]
    # output_files = ["results/%s/cnv_sv/%s.svdb_merged.vcf.gz" % (sample, sample) for sample in get_samples(samples)]
    # output_files = ["results/%s/%s.contamination.html" % (sample, sample) for sample in get_samples(samples)]
    # output_files = ["results/%s/expansionhunter_reviewer/" % (sample) for sample in get_samples(samples)]
    # output_files = ["results/%s/%s.expansionhunter_stranger.vcf.gz" % (sample, sample) for sample in get_samples(samples)]
    # output_files = ["results/%s/%s.coverage_analysis.xlsx" % (sample, sample) for sample in get_samples(samples)]
    # output_files = ["results/%s/SMNCopyNumberCaller/%s.smn_charts.pdf" % (sample, sample) for sample in get_samples(samples)]
    # output_files = ["results/%s/SMNCopyNumberCaller/%s.smn_caller.tsv" % (sample, sample) for sample in get_samples(samples)]
    # output_files = ["results/%s/SMNCopyNumberCaller/%s.smn_caller.json" % (sample, sample) for sample in get_samples(samples)]
    # output_files = [
    #     "results/%s/%s_%s.crumble.cram.crai" % (sample, sample, unit_type)
    #     for sample in get_samples(samples)
    #     for unit_type in get_unit_types(units, sample)
    # ]
    # output_files = [
    #     "results/%s/%s_%s.crumble.cram" % (sample, sample, unit_type)
    #     for sample in get_samples(samples)
    #     for unit_type in get_unit_types(units, sample)
    # ]
    # output_files = ["results/%s/%s_snv_indels.vcf.gz.tbi" % (sample, sample) for sample in get_samples(samples)]
    # output_files = ["results/%s/%s_snv_indels.vcf.gz" % (sample, sample) for sample in get_samples(samples)]
    # output_files = ["results/%s/%s_snv_indels.filtered.vcf.gz.tbi" % (sample, sample) for sample in get_samples(samples)]
    # output_files = ["results/%s/%s_snv_indels.filtered.vcf.gz" % (sample, sample) for sample in get_samples(samples)]
    return output_files
### Include copy all files we want to transfer
# def compile_output_list(wildcards):
#     output_files = []
#     types = set([unit.type for unit in units.itertuples()])
#     for output in output_json:
#         output_files = set(
#             [
#                 output.format(sample=sample, type=unit_type, caller=caller)
#                 for sample in get_samples(samples)
#                 for unit_type in get_unit_types(units, sample)
#                 if unit_type in set(output_json[output]["types"]).intersection(types)
#                 for caller in config["bcbio_variation_recall_ensemble"]["callers"]
#             ]
#         )
#     return list(set(output_files))
#
# def generate_copy_code(workflow, output_json):
#     code = ""
#     for result, values in output_json.items():
#         if values["file"] is not None:
#             input_file = values["file"]
#             output_file = result
#             rule_name = values["name"]
#             mem_mb = config.get('_copy', {}).get("mem_mb", config["default_resources"]["mem_mb"])
#             mem_per_cpu = config.get('_copy', {}).get("mem_mb", config["default_resources"]["mem_mb"])
#             partition = config.get("_copy", {}).get("partition", config["default_resources"]["partition"])
#             threads = config.get("_copy", {}).get("threads", config["default_resources"]["threads"])
#             time = config.get("_copy", {}).get("time", config["default_resources"]["time"])
#             copy_container = config.get("_copy", {}).get("container", config["default_container"])
#             result_file = os.path.basename(output_file)
#             code = f'@workflow.rule(name="{rule_name}")\n'
#             code = f'@workflow.input("{input_file}")\n'
#             code = f'@workflow.output("{output_file}")\n'
#             code = f'@workflow.log("logs/{rule_name}_{result_file}.log")\n'
#             code = f'@workflow.container("{copy_container}")\n'
#             code = f'@workflow.conda("../env/copy_result.yaml")\n'
#             code = f'@workflow.resources(time = "{time}", threads = {threads}, mem_mb = {mem_mb}, mem_per_cpu = {mem_per_cpu}, partition = "{partition}")\n'
#             code = '@workflow.shellcmd("cp {input} {output}")\n\n'
#             code = "@workflow.run\n"
#             code = (
#                 f"def __rule_{rule_name}(input, output, params, wildcards, threads, resources, log, version, rule, "
#                 "conda_env, container_img, singularity_args, use_singularity, env_modules, bench_record, jobid, is_shell, "
#                 "bench_iteration, cleanup_scripts, shadow_dir, edit_notebook, conda_base_path, basedir, runtime_sourcecache_path, "
#                 "__is_snakemake_rule_func=True):\n"
#                 '\tshell ( "(cp {input[0]} {output[0]}) &> {log}" , bench_record=bench_record, bench_iteration=bench_iteration)\n\n'
#             )
#     exec(compile(code, "result_to_copy", "exec"), workflow.globals)
#
#
# generate_copy_code(workflow, output_json)

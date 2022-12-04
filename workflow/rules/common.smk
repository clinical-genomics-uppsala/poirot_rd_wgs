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

### Read and validate units file
units = (
    pandas.read_table(config["units"], dtype=str)
    .set_index(["sample", "type", "flowcell", "lane", "barcode"], drop=False)
    .sort_index()
)
validate(units, schema="../schemas/units.schema.yaml")


### Set wildcard constraints
wildcard_constraints:
    barcode="[A-Z+]+",
    chr="[^_]+",
    flowcell="[A-Z0-9]+",
    lane="L[0-9]+",
    read="fastq[1|2]",
    sample="|".join(get_samples(samples)),
    type="N|T|R",


### Functions
def get_flowcell(units, wildcards):
    flowcells = set([u.flowcell for u in get_units(units, wildcards)])
    if len(flowcells) > 1:
        raise ValueError("Sample type combination from different sequence flowcells")
    return flowcells.pop()


def get_in_fastq(units, wildcards):
    return expand(
        "prealignment/merged/{{sample}}_{{type}}_{read}.fastq.gz",
        read=["fastq1", "fastq2"],
    )


def get_input_fastq(units, wildcards):
    return expand(
        "prealignment/fastp_pe/{{sample}}_{flowcell_lane_barcode}_{{type}}_{read}.fastq.gz",
        flowcell_lane_barcode=[
            "{}_{}_{}".format(unit.flowcell, unit.lane, unit.barcode) for unit in get_units(units, wildcards, wildcards.type)
        ],
        read=["fastq1", "fastq2"],
    )


def get_in_fq(wildcards):
    input_list = []
    for unit in get_units(units, wildcards, wildcards.type):
        prefix = "prealignment/fastp_pe/{}_{}_{}_{}_{}".format(unit.sample, unit.flowcell, unit.lane, unit.barcode, unit.type)
        input_unit = "{}_fastq1.fastq.gz {}_fastq2.fastq.gz {}".format(
            prefix,
            prefix,
            "'@RG\\tID:{}\\tSM:{}\\tPL:{}\\tPU:{}\\tLB:{}'".format(
                "{}_{}.{}.{}".format(unit.sample, unit.type, unit.lane, unit.barcode),
                "{}_{}".format(unit.sample, unit.type),
                unit.platform,
                "{}.{}.{}".format(unit.flowcell, unit.lane, unit.barcode),
                "{}_{}".format(unit.sample, unit.type),
            ),
        )
        input_list.append(input_unit)
    return " --in-fq ".join(input_list)


def get_in_gvcf(wildcards):
    gvcf_list = [
        "snv_indels/deepvariant/{}_{}.g.vcf".format(sample, t)
        for sample in get_samples(samples)
        for t in get_unit_types(units, sample)
    ]
    return " -i ".join(gvcf_list)


def get_peddy_sex(wildcards, peddy_sex_check):
    sample = '{}_{}'.format(wildcards.sample, wildcards.type)
    sex_df = pd.read_table(peddy_sex_check, sep=',').set_index("sample_id", drop=False)

    sample_sex = sex_df.at[sample, 'predicted_sex']

    return sample_sex


def get_locus_str(loci):
    with open(loci, 'r') as catfile:
        loc_str = catfile.readline().rstrip()
    return loc_str


def compile_output_list(wildcards: snakemake.io.Wildcards):
    files = {
        "cnv_sv/cnvpytor": ["vcf"],
        "cnv_sv/expansionhunter": ["vcf"],
        "cnv_sv/stranger": ["stranger.vcf"],
        "cnv_sv/svdb_query": ["svdb_query.vcf"],
        "cnv_sv/tiddit": ["vcf"],
        "compression/crumble": ["crumble.cram"],
        "qc/create_cov_excel": ["coverage.xlsx"],
    }
    output_files = [
        "%s/%s_%s.%s" % (prefix, sample, unit_type, suffix)
        for prefix in files.keys()
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        for suffix in files[prefix]
    ]
    output_files += [
        "cnv_sv/expansionhunter/reviewer/%s_%s/" % (sample, unit_type)
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
    ]
    output_files += [
        "cnv_sv/manta_run_workflow_n/%s/results/variants/diploidSV.vcf.gz" % (sample) for sample in get_samples(samples)
    ]
    output_files += ["qc/multiqc/multiqc_DNA.html"]
    output_files += [
        "qc/peddy/peddy.peddy.ped",
        "qc/peddy/peddy.ped_check.csv",
        "qc/peddy/peddy.sex_check.csv",
        "qc/peddy/peddy.het_check.csv",
        "qc/peddy/peddy.html",
        "qc/peddy/peddy.vs.html",
        "qc/peddy/peddy.background_pca.json",
    ]
    output_files += [
        "compression/spring/%s_%s_%s_%s_%s.spring" % (sample, flowcell, lane, barcode, t)
        for sample in set(units["sample"])
        for flowcell in set(units["flowcell"])
        for lane in set(units["lane"])
        for barcode in set(units["barcode"])
        for t in set(units["type"])
    ]
    output_files += ["vcf_final/%s.vcf.gz.tbi" % (sample) for sample in get_samples(samples)]
    return output_files



### Include copy all files we want to transfer
# def compile_output_list(wildcards):
#     output_files = []
#     types = set([unit.type for unit in units.itertuples()])
#     for output in output_json:
#         output_files += set(
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
#             code += f'@workflow.rule(name="{rule_name}")\n'
#             code += f'@workflow.input("{input_file}")\n'
#             code += f'@workflow.output("{output_file}")\n'
#             code += f'@workflow.log("logs/{rule_name}_{result_file}.log")\n'
#             code += f'@workflow.container("{copy_container}")\n'
#             code += f'@workflow.conda("../env/copy_result.yaml")\n'
#             code += f'@workflow.resources(time = "{time}", threads = {threads}, mem_mb = {mem_mb}, mem_per_cpu = {mem_per_cpu}, partition = "{partition}")\n'
#             code += '@workflow.shellcmd("cp {input} {output}")\n\n'
#             code += "@workflow.run\n"
#             code += (
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

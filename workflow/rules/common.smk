___author__ = "Jessika Nordin, Padraic Corcoran"
__copyright__ = "Copyright 2022"
__email__ = "jessika.nordin@scilifelab.uu.se"
__license__ = "GPL-3"

import pandas
import yaml
from datetime import datetime

from hydra_genetics.utils.misc import get_module_snakefile
from hydra_genetics.utils.resources import load_resources
from hydra_genetics.utils.misc import replace_dict_variables
from hydra_genetics.utils.samples import *
from hydra_genetics.utils.units import *
from hydra_genetics.utils.misc import extract_chr
from snakemake.utils import min_version
from snakemake.utils import validate
from hydra_genetics import min_version as hydra_min_version

from hydra_genetics.utils.misc import export_config_as_file
from hydra_genetics.utils.software_versions import add_version_files_to_multiqc
from hydra_genetics.utils.software_versions import add_software_version_to_config
from hydra_genetics.utils.software_versions import export_pipeline_version_as_file
from hydra_genetics.utils.software_versions import export_software_version_as_file
from hydra_genetics.utils.software_versions import get_pipeline_version
from hydra_genetics.utils.software_versions import touch_pipeline_version_file_name
from hydra_genetics.utils.software_versions import touch_software_version_file
from hydra_genetics.utils.software_versions import use_container


hydra_min_version("3.0.0")

min_version("7.8.0")

if not workflow.overwrite_configfiles:
    sys.exit("At least one config file must be passed using --configfile/--configfiles, by command line or a profile!")

try:
    config = replace_dict_variables(config)
    validate(config, schema="../schemas/config.schema.yaml")
except WorkflowError as we:
    # Probably a validation error, but the original exsception in lost in
    # snakemake. Pull out the most relevant information instead of a potentially
    # *very* long error message.
    if not we.args[0].lower().startswith("error validating config file"):
        raise
    error_msg = "\n".join(we.args[0].splitlines()[:2])
    parent_rule_ = we.args[0].splitlines()[3].split()[-1]
    if parent_rule_ == "schema:":
        sys.exit(error_msg)
    else:
        schema_hiearachy = parent_rule_.split()[-1]
        schema_section = ".".join(re.findall(r"\['([^']+)'\]", schema_hiearachy)[1::2])
        sys.exit(f"{error_msg} in {schema_section}")

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

## read the output json
with open(config["output"]) as output:
    output_json = json.load(output)

## get version information on pipeline, containers and software

pipeline_name = "Poirot"
pipeline_version = get_pipeline_version(workflow, pipeline_name=pipeline_name)
version_files = touch_pipeline_version_file_name(
    pipeline_version, date_string=pipeline_name, directory="results/versions/software"
)
if use_container(workflow):
    version_files.append(touch_software_version_file(config, date_string=pipeline_name, directory="results/versions/software"))
add_version_files_to_multiqc(config, version_files)


onstart:
    export_pipeline_version_as_file(pipeline_version, date_string=pipeline_name, directory="results/versions/software")
    if use_container(workflow):
        update_config, software_info = add_software_version_to_config(config, workflow, False)
        export_software_version_as_file(software_info, date_string=pipeline_name, directory="results/versions/software")
    date_string = datetime.now().strftime("%Y%m%d")
    export_config_as_file(update_config, date_string=date_string, directory="results/versions")


### Set wildcard constraints
wildcard_constraints:
    barcode="[A-Z+]+",
    chr="[^_]+",
    flowcell="[A-Z0-9]+",
    lane="L[0-9]+",
    read="fastq[1|2]",
    sample="|".join(get_samples(samples)),
    type="N|T|R",
    vcf="vcf|g.vcf|unfiltered.vcf",


## contigs in hg38
contigs = extract_chr("%s.fai" % (config.get("reference", {}).get("fasta", "")), filter_out=[])
skip_contigs = [c for c in contigs if "_" in c or c == "chrEBV"]

### Functions


def get_bam_input(wildcards, use_sample_wildcard=True, use_type_wildcard=True):
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
        # if by_chr:  # if a bam for single chromosome is needed
        #     bam_input = "alignment/picard_mark_duplicates/{}_{}.bam".format(sample_str, wildcards.chr)
        bam_input = "alignment/samtools_merge_bam/{}.bam".format(sample_str)
    else:
        sys.exit("valid options for aligner are: bwa_gpu or bwa_cpu")

    bai_input = "{}.bai".format(bam_input)

    return (bam_input, bai_input)


def get_chrom_deepvariant_vcfs(wildcards, vcf_type):
    skip_contig_patterns = config.get("reference", {}).get("merge_contigs", "")
    skip_contigs = []
    ref_fasta = config.get("reference", {}).get("fasta", "")
    all_contigs = extract_chr(f"{ref_fasta}.fai", filter_out=[])
    for pattern in contig_patterns:
        for contig in all_contigs:
            contig_match = re.match(pattern, contig)
            if contig_match is not None:
                skip_contigs.append(contig_match.group())

    chroms = extract_chr(f"{ref_fasta}.fai", filter_out=skip_contigs)
    vcf_suffix = "vcf.gz"
    if vcf_type == "gvcf":
        vcf_suffix = "g.vcf.gz"

    vcf_list = [f"snv_indels/deepvariant/{wildcards.sample}_{wildcards.type}_{chr}.{vcf_suffix}" for chr in chroms]
    tbi_list = [f"{v}.tbi" for v in vcf_list]

    return (vcf_list, tbi_list)


def get_vcf_input(wildcards):
    caller = config.get("snp_caller", None)
    if caller is None:
        sys.exit("snp_caller missing from config, valid options: deepvariant_gpu or deepvariant_cpu")
    elif caller == "deepvariant_gpu":
        vcf_input = f"parabricks/pbrun_deepvariant/{wildcards.sample}_{wildcards.type}.vcf"
    elif caller == "deepvariant_cpu":
        vcf_input = f"snv_indels/deepvariant/{wildcards.sample}_{wildcards.type}.merged.vcf"
    else:
        sys.exit("Invalid options for snp_caller, valid options are: deepvariant_gpu or deepvariant_cpu")

    return vcf_input


def get_gvcf_list(wildcards):
    caller = config.get("snp_caller", None)
    if caller is None:
        sys.exit("snp_caller missing from config, valid options: deepvariant_gpu or deepvariant_cpu")
    elif caller == "deepvariant_gpu":
        gvcf_path = "parabricks/pbrun_deepvariant"
        gvcf_list = [
            "{}/{}_{}.g.vcf".format(gvcf_path, sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    elif caller == "deepvariant_cpu":
        gvcf_path = "snv_indels/deepvariant"
        gvcf_list = [
            "{}/{}_{}.merged.g.vcf".format(gvcf_path, sample, t)
            for sample in get_samples(samples)
            for t in get_unit_types(units, sample)
        ]
    else:
        sys.exit("Invalid options for snp_caller, valid options are: deepvariant_gpu or deepvariant_cpu")

    return gvcf_list


def get_gvcf_trio(wildcards):
    caller = config.get("snp_caller", None)

    proband_sample = samples[samples.index == wildcards.sample]
    trio_id = proband_sample.at[wildcards.sample, "trioid"]
    mother_sample = samples[(samples.trio_member == "mother") & (samples.trioid == trio_id)].index[0]
    father_sample = samples[(samples.trio_member == "father") & (samples.trioid == trio_id)].index[0]

    if caller is None:
        sys.exit("snp_caller missing from config, valid options: deepvariant_gpu or deepvariant_cpu")
    elif caller == "deepvariant_gpu":
        child_gvcf = "parabricks/pbrun_deepvariant/{}_{}.g.vcf".format(wildcards.sample, wildcards.type)
        mother_gvcf = "parabricks/pbrun_deepvariant/{}_{}.g.vcf".format(mother_sample, wildcards.type)
        father_gvcf = "parabricks/pbrun_deepvariant/{}_{}.g.vcf".format(father_sample, wildcards.type)
    elif caller == "deepvariant_cpu":
        child_gvcf = "snv_indels/deepvariant/{}_{}.merged.g.vcf".format(wildcards.sample, wildcards.type)
        mother_gvcf = "snv_indels/deepvariant/{}_{}.merged.g.vcf".format(mother_sample, wildcards.type)
        father_gvcf = "snv_indels/deepvariant/{}_{}.merged.g.vcf".format(father_sample, wildcards.type)
    else:
        sys.exit("Invalid options for snp_caller, valid options are: deepvariant_gpu or deepvariant_cpu")

    gvcf_list = [child_gvcf, mother_gvcf, father_gvcf]

    return gvcf_list


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
    gvcf_input = "-i {}".format(" -i ".join(input.gvcfs))

    return gvcf_input


def get_vcfs_for_svdb_merge(wildcards, input):
    vcfs_with_suffix = []
    vcfs_with_suffix.append(f"{input.manta}:manta")
    vcfs_with_suffix.append(f"{input.cnvpytor}:cnvpytor")

    return vcfs_with_suffix


def get_str_panel_list(wildcards):

    panel_dir = config["reference"]["str_panels_dir"]
    panel = wildcards.panel
    panel_list_path = f"{panel_dir}/{panel}.list"

    return panel_list_path


def compile_output_list(wildcards):
    output_files = []
    types = set([unit.type for unit in units.itertuples()])
    for output in output_json:
        if output == "results/{sample}/{sample}.upd_regions.bed":
            for sample in samples[samples.trio_member == "proband"].index:
                proband_trio_id = samples[samples.index == sample].trioid.iloc[0]
                trio_num = samples[samples.trioid == proband_trio_id].shape[0]
                if trio_num != 3:
                    continue
                else:
                    output_files.append(output.format(sample=sample))
        else:
            output_files += set(
                [
                    output.format(sample=sample, flowcell=flowcell, lane=lane, barcode=barcode, type=unit_type, panel=str_panel)
                    for sample in get_samples(samples)
                    for unit_type in get_unit_types(units, sample)
                    if unit_type in set(output_json[output]["types"])
                    for flowcell in set([u.flowcell for u in units.loc[(sample, unit_type)].dropna().itertuples()])
                    for barcode in set([u.barcode for u in units.loc[(sample, unit_type)].dropna().itertuples()])
                    for lane in set([u.lane for u in units.loc[(sample, unit_type)].dropna().itertuples()])
                    for str_panel in [panel_list.split(".")[0] for panel_list in config["reference"]["str_panels"]]
                ]
            )

    return list(set(output_files))


def generate_copy_code(workflow, output_json):
    code = ""
    for result, values in output_json.items():
        if values["file"] is not None:
            input_file = values["file"]
            output_file = result
            rule_name = values["name"]
            mem_mb = config.get("_copy", {}).get("mem_mb", config["default_resources"]["mem_mb"])
            mem_per_cpu = config.get("_copy", {}).get("mem_mb", config["default_resources"]["mem_mb"])
            partition = config.get("_copy", {}).get("partition", config["default_resources"]["partition"])
            threads = config.get("_copy", {}).get("threads", config["default_resources"]["threads"])
            time = config.get("_copy", {}).get("time", config["default_resources"]["time"])
            copy_container = config.get("_copy", {}).get("container", config["default_container"])
            code += f'@workflow.rule(name="{rule_name}")\n'
            code += f'@workflow.input("{input_file}")\n'
            if rule_name == "_copy_reviewer":  # handle rule that has directory as output
                result_file = "{sample}"
                code += f'@workflow.output(directory("{output_file}"))\n'
            else:
                result_file = os.path.basename(output_file)
                code += f'@workflow.output("{output_file}")\n'
            code += f'@workflow.log("logs/{rule_name}_{result_file}.log")\n'
            code += f'@workflow.container("{copy_container}")\n'
            code += f'@workflow.resources(time = "{time}", threads = {threads}, mem_mb = {mem_mb}, mem_per_cpu = {mem_per_cpu}, partition = "{partition}")\n'
            code += '@workflow.shellcmd("cp -r {input} {output}")\n\n'
            code += "@workflow.run\n"
            code += (
                f"def __rule_{rule_name}(input, output, params, wildcards, threads, resources, log, version, rule, "
                "conda_env, container_img, singularity_args, use_singularity, env_modules, bench_record, jobid, is_shell, "
                "bench_iteration, cleanup_scripts, shadow_dir, edit_notebook, conda_base_path, basedir, runtime_sourcecache_path, "
                "__is_snakemake_rule_func=True):\n"
                '\tshell ( "(cp -r {input[0]} {output[0]}) &> {log}" , bench_record=bench_record, bench_iteration=bench_iteration)\n\n'
            )

    exec(compile(code, "result_to_copy", "exec"), workflow.globals)


generate_copy_code(workflow, output_json)

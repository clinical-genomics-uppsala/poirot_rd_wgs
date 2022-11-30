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
        "cnv_sv/tiddit": ["vcf"],
        "cnv_sv/svdb_query": ["svdb_query.vcf"],
        "mitochondrial/gatk_split_multi_allelic_sites": ["vcf"]
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
    output_files += ["vcf_final/%s.vcf.gz.tbi" % (sample) for sample in get_samples(samples)]
    return output_files

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
    flowcell="[A-Z0-9]+",
    lane="L[0-9]+",
    sample="|".join(get_samples(samples)),
    unit="N|T|R",
    read="fastq[1|2]",


### Functions
def get_in_fastq(units, wildcards):
    return expand(
        "prealignment/merged/{{sample}}_{{type}}_{read}.fastq.gz",
        read=["fastq1", "fastq2"],
    )


def get_in_fq(wildcards):
    input_list = []
    for unit in get_units(units, wildcards, wildcards.type):
        prefix = "prealignment/merged/{}_{}".format(unit.sample, unit.type)
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


def compile_output_list(wildcards: snakemake.io.Wildcards):
    files = {
        #        "cnv_sv/cnvpytor": ["vcf"],
#        "cnv_sv/expansionhunter": ["vcf", "stranger.vcf"],
        "cnv_sv/tiddit": ["vcf"],
#        "vcf_final": ["vcf.gz.tbi"],
    }

    output_files = [
        "%s/%s_%s.%s" % (prefix, sample, unit_type, suffix)
        for prefix in files.keys()
        for sample in get_samples(samples)
        for unit_type in get_unit_types(units, sample)
        for suffix in files[prefix]
    ]
#    output_files += [
#        "cnv_sv/expansionhunter/reviewer/%s_%s/" % (sample, unit_type)
#        for sample in get_samples(samples)
#        for unit_type in get_unit_types(units, sample)
#    ]
    output_files += ["cnv_sv/manta_run_workflow_n/%s/results/variants/candidateSV.vcf.gz" % (sample)
        for sample in get_samples(samples)
    ]
#    output_files += ["qc/multiqc/multiqc_DNA.html"]
    return output_files

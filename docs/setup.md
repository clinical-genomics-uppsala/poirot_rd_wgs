# Pipeline setup and configuration
There are a number of main files that governs how the pipeline is executed listed below:

* Snakefile
* common.smk
* config.yaml
* resources.yaml
* profile/uppsala/config.yaml
* samples.tsv and units.tsv

There is more general information about the content of these files in hydra-genetics documentation in [code standards](https://hydra-genetics.readthedocs.io/en/latest/development/standards/), [config](https://hydra-genetics.readthedocs.io/en/latest/development/standards/#config) and [Snakefile](https://hydra-genetics.readthedocs.io/en/latest/make_pipeline/import/).

## Snakefile
The `Snakefile` is located in workflow/ and imports hydra-genetics modules and rules as well as modifies these rules when needed. It also imports pipeline specific rules and define rule orders. Finally, this is where the rule all is defined.

## common.smk
The `common.smk` is located under workflow/rules/. This is a general rule taking care of any actions that are not directly connected with running a specific program. It includes version checks, import of config, resources, tsv-files and validations using schemas. Functions used by pipeline specific rules are also defined here as well as the output files using the function **compile_output_list** which programmatically generates a list of all necessary output files for the module to be targeted in the all rule defined in the `Snakemake` file. See further [Result files](https://hydra-genetics.readthedocs.io/en/latest/make_pipeline/results/).

## config.yaml
The `config.yaml` is located under config/. The file ties all file and other dependencies as well as parameters for different rules together.
See further [pipeline configuration](https://hydra-genetics.readthedocs.io/en/latest/make_pipeline/config/).

<br />

/// details | Expand to view current config.yaml
```yaml
{% include "includes/config.yaml" %}
```
///


## resources.yaml
The `resources.yaml` is located under config/. The file declares default resources used by rules as well as resources for specific rules that needs more resources than allocated by default. See further [pipeline configuration](https://hydra-genetics.readthedocs.io/en/latest/make_pipeline/config/).

```yaml
# ex, default resources
default_resources:
  threads: 1
  time: "4:00:00"
  mem_mb: 6144
  mem_per_cpu: 6144
  partition: "low"

# ex, rule override
vardict:
  time: "48:00:00"
```

<br />

/// details | Expand to view current resources.yaml
```yaml
{% include "includes/resources.yaml" %}
```
///

## profile yaml
Profiles are saved in yaml files and used to control how snakemake will be executed, if jobs will be submitted
to a cluster, use singularity, restart on failure and so forth. It also forward requested resources to drmaa using
a drmaa variable.

```yaml
# ex, snakemake settings
jobs: 100
keep-going: True
restart-times: 2
rerun-incomplete: True
use-singularity: True
configfile: "config/config.yaml"
singularity-args: "-e --cleanenv -B /projects -B /data -B /beegfs
```

```yaml
# ex, drmaa settings
drmaa: " -A wp1 -N 1-1 -t {resources.time} -n {resources.threads} --mem={resources.mem_mb} --mem-per-cpu={resources.mem_per_cpu} --mem-per-cpu={resources.mem_per_cpu} --partition={resources.partition} -J {rule} -e slurm_out/{rule}_%j.err -o slurm_out/{rule}_%j.out"
drmaa-log-dir: "slurm_out"
default-resources: [threads=1, time="04:00:00", partition="low", mem_mb="3074", mem_per_cpu="3074"]
```

## samples.tsv and units.tsv
The `samples.tsv` and `units.tsv` are input files that must be generated before running the pipeline and should in general be located in the base folder of the analysis folder, can be changed in the config.yaml. See further [running the pipeline](running.md) and [create input files](https://hydra-genetics.readthedocs.io/en/latest/run_pipeline/create_sample_files/).

### Example samples.tsv

| sample | tumor_content | sex	| trioid	| trio_member	|
|-|-| -	| -	| -	|
| NA12878	| 0.0	| female	| NA	| NA	|
| NA12911	| 0	| male	| NA	| NA	|

### Example units.tsv

| sample | type | machine | platform | flowcell | lane | barcode | fastq1 | fastq2 | adapter |
|-|-|-|-|-|-|-|-|-|-|
| NA12878	| N	| NDX550407_RUO	| NextSeq	| HKTG2BGXG	| L001 | ACGGAACA+ACGAGAAC | fastq/NA12878_fastq1.fastq.gz | fastq/NA12878_fastq2.fastq.gz | ACGT,ACGT |
| NA12911	| N	| NDX550407_RUO	| NextSeq	| HKTG2BGXG	| L001 | TCGGAACT+TCGAGAAT | fastq/NA12911_fastq1.fastq.gz | fastq/NA12911_fastq2.fastq.gz | ACGT,ACGT |


## Coverage analysis
Too get coverage analysis for gene panels three things are needed 
  1) a file with a list of gene panels and 
  2) a list of which genes that is in the panel and 
  3) a bed file with information about the genes.

### Gene panels
There are several gene panels available at [Poirot's config git](https://github.com/clinical-genomics-uppsala/poirot_config/tree/main/config/gene_panels).

#### 1) example genepanel.list
BRCA

CADASIL

EBS

EDS

#### 2) example gene panel list BRCA.list
APC

ATM

BAP1

BMPR1A

BRCA1

BRCA2

### 3) Bed-file with gene information

#### Example bed-file

| chr1 | 999056   | 999434 	| HES4_NM_021170.4_[4]	|

| chr1 | 999523   | 999615	| HES4_NM_021170.4_[3]	|

| chr1 | 999689   | 999789	| HES4_NM_021170.4_[2]	|

| chr1 | 999863   | 999975	| HES4_NM_021170.4_[1]	|

| chr1 | 1013571  | 1013578	| ISG15_NM_005101.4_[1]	|

| chr1 | 1013981  | 1014480	| ISG15_NM_005101.4_[2]	|

<br />


<br />

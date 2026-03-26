# Rules specific to Poirot

## create_cov_excel
Script that creates the gene panel coverage excel file for poirot.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__coverage__create_cov_excel#

#### :left_right_arrow: input / output files


---

## deepvariant_add_ref
Add the reference genome path to the deepvariant vcf header

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__add_ref_to_vcf__deepvariant_add_ref#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__add_ref_to_vcf__deepvariant_add_ref#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__deepvariant_add_ref#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__deepvariant_add_ref#

--- 
## [filter_par_dups]
A custom python script to filter DUP calls in male sample
chrX PAR regions in cnvpytor vcf files.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__filter_par_dups__filter_par_dups#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__filter_par_dups__filter_par_dups#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__filter_par_dups#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__filter_par_dups#

---

## fix_mt_gt
A script that postprocesses the GATK mitochondrial normalised mutect2 VCF. It looks for GT fields that have more than two entries (e.g. '0/././1, or '0/1/./.' etc)
and converts them to '0/1' as some tools can not parse the vcf when the 
GT field has missing alleles and has more than two allele fields.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__fix_mt_gt__fix_mt_gt#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__fix_mt_gt__fix_mt_gt#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__fix_mt_gt#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__fix_mt_gt#

---

## svdb_add_ref
Add the reference genome path to svdb merge vcf header

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__add_ref_to_vcf__svdb_add_ref#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__add_ref_to_vcf__svdb_add_ref#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__svdb_add_ref#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__svdb_add_ref#

---

## tiddit_add_ref
Add the reference genome path to tiddit vcf header

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__add_ref_to_vcf__tiddit_add_ref#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__add_ref_to_vcf__tiddit_add_ref#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__tiddit_add_ref#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__tiddit_add_ref#

---

## vcf_to_aed
Conversion of cnvpytor vcf to AED file format. The AED file can be read by Chromosome Analysis Suite (ChAS).

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__vcf_to_aed__vcf_to_aed#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__vcf_to_aed__vcf_to_aed#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__vcf_to_aed#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__vcf_to_aed#

---

## vcf_to_aed_filtered
Conversion of teh filtered cnvpytor vcf to AED file format. The AED file can be read by Chromosome Analysis Suite (ChAS).

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__vcf_to_aed__vcf_to_aed_filtered#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__vcf_to_aed__vcf_to_aed_filtered#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__vcf_to_aed_filtered#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__vcf_to_aed_filtered#

---

## create_somalier_mqc_tsv
Create MultiQC custom content TSV files from Somalier output. This script processes somalier relatedness and sex check data to create custom tables similar to Peddy tables, with Pass/Fail QC checks.

### :snake: Rule

#SNAKEMAKE_RULE_SOURCE__somalier__create_somalier_mqc_tsv#

#### :left_right_arrow: input / output files

#SNAKEMAKE_RULE_TABLE__somalier__create_somalier_mqc_tsv#

### :wrench: Configuration

#### Software settings (`config.yaml`)

#CONFIGSCHEMA__create_somalier_mqc_tsv#

#### Resources settings (`resources.yaml`)

#RESOURCESSCHEMA__create_somalier_mqc_tsv#



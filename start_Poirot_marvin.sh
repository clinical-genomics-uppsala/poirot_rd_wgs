#!/usr/bin/env bash

runfolder=$1
runname=$(basename $1)
year=$(basename $(dirname $1))
analysis_type=$(basename $(dirname $(dirname $(dirname $1))))

echo $runfolder
echo $runname
echo $year
echo $analysis_type

module add slurm-drmaa
module add singularity
module add snakemake

cd ${runfolder};

mv SampleSheet.csv ${runname}_samplesheet.csv;

snakemake --profile marvin -s workflow/Snakefile --configfile config/config.yaml


#### Not done yet!!! ####

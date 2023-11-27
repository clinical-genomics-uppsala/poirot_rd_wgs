#!/bin/bash

set -e
set -u
set -o pipefail

output=/beegfs-storage/projects/wp3/nobackup/WGS/NovaSeqX_validation_runs/20231030_LH00338_0003_A22FJ27LT3/poirot_output

if [ ! -d ${output} ]; then
    mkdir ${output}
    mkdir ${output}/slurm
fi

# source virtual env with requirements.txt installed
python3.9 -m venv hydra_env
source hydra_env/bin/activate
pip install -r requirements.txt 

module load singularity/3.7.1
module load slurm-drmaa/1.1.3 

sample_sheet=/beegfs-storage/projects/wp3/nobackup/WGS/NovaSeqX_validation_runs/20231030_LH00338_0003_A22FJ27LT3/SampleSheet.csv
fastq_dir=/beegfs-storage/projects/wp3/nobackup/WGS/NovaSeqX_validation_runs/20231030_LH00338_0003_A22FJ27LT3/fastq/

## barcodes not currently in the Nextseq fastq headers, so use Ns 
hydra-genetics create-input-files -d ${fastq_dir} -t N --tc 0  -f

# adds trio info to samples.tsv and renames to SampleSheet.csv ids, create sample order files for multiqc
python extract_sample_sheet_info.py ${sample_sheet} 

rsync -ruv config ${output}

snakemake  --profile profiles/slurm/ --configfile config/config.yaml -p  -d ${output} --notemp

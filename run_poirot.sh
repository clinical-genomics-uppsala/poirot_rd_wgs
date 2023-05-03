#!/bin/bash

set -e
set -u
set -o pipefail

# source virtual env with requirements.txt installed
source ~/hydra-env/bin/activate 

module load singularity/3.7.1
module load slurm-drmaa/1.1.3 

sample_sheet=/beegfs-storage/projects/wp3/nobackup/WGS/WGS_pilot_CGU202111/THG11/230330_A00181_0643_AHNN37DSX5/SampleSheet.csv
fastq_dir=/beegfs-storage/projects/wp3/nobackup/WGS/WGS_pilot_CGU202111/THG11/230330_A00181_0643_AHNN37DSX5/

## barcodes not currently in the Nextseq fastq headers, so use Ns 
hydra-genetics create-input-files -d ${fastq_dir} -t N --tc 0  -f

# adds trio info to samples.tsv and renames to SampleSheet.csv ids, create sample order files for multiqc
python extract_sample_sheet_info.py ${sample_sheet} 

snakemake  --profile profiles/slurm/ --configfile config/config.yaml -pn

#!/bin/bash

set -e
set -u
set -o pipefail

# source virtual env with requirements.txt installed
python3.9 -m venv hydra_env
source hydra_env/bin/activate
pip install -r requirements.txt 

module load singularity/3.7.1
module load slurm-drmaa/1.1.3 

sample_sheet=/beegfs-storage/projects/wp4/nobackup/DDS_downloads/DataDelivery_2023-08-22_16-10-05_snpseq00449/files/WD-3649/230818_A00181_0673_AHKCL2DSX7/SampleSheet.csv
fastq_dir=/beegfs-storage/projects/wp4/nobackup/DDS_downloads/DataDelivery_2023-08-22_16-10-05_snpseq00449/files/WD-3649/230818_A00181_0673_AHKCL2DSX7/

## barcodes not currently in the Nextseq fastq headers, so use Ns 
hydra-genetics create-input-files -d ${fastq_dir} -t N --tc 0  -f

# adds trio info to samples.tsv and renames to SampleSheet.csv ids, create sample order files for multiqc
python extract_sample_sheet_info.py ${sample_sheet} 

snakemake  --profile profiles/slurm/ --configfile config/config.yaml -p

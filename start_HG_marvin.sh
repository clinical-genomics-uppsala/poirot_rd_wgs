#!/usr/bin/env bash
# To run script:
# bash /projects/wp3/nobackup/TWIST/Bin/Poirot/start_HG_marvin.sh /projects/wp3/nobackup/WGS/*/230619_NDX550925_RUO_0133_AHJ2YMBGXT/Sample_* HG4

set -euo pipefail

module load slurm-drmaa/1.1.3
module load singularity/3.7.1

poirotFolder=/beegfs-storage/projects/wp3/nobackup/TWIST/Bin/Poirot

python3.9 -m venv ${poirotFolder}/hydra_env
source ${poirotFolder}/hydra_env/bin/activate
pip install -r ${poirotFolder}/requirements.txt

fastqFolder=$1
sequencerun=$2    #Sequence ID
startDir=$PWD

outbox=$(echo $(echo ${poirotFolder} | rev | cut -d/ -f3- | rev)/OUTBOX)

# Create outbox and scratch folders
echo 'Creating outbox and scratch folders'
if [ ! -d "/scratch/wp3/TWIST/${sequencerun}/" ]
then
  mkdir /scratch/wp3/TWIST/${sequencerun}/
fi

if [ ! -d "${outbox}/${sequencerun}/" ]
then
  mkdir ${outbox}/${sequencerun}/
fi

# Cp data to scratch
echo 'Copy to scratch' && \
#To save space on scratch I comment out copy fastq-files there
#rsync -ru ${fastqFolder%}/Sample_WD* /scratch/wp3/TWIST/${sequencerun}/fastq  && \
rsync -ru SampleSheet.csv /scratch/wp3/TWIST/${sequencerun}/  && \
rsync -ru ${poirotFolder}/config /scratch/wp3/TWIST/${sequencerun}/  && \

cd /scratch/wp3/TWIST/${sequencerun}/  && \

# Create sample and unit files
#To save space on scratch we leave the fastq-files in nobackup
#hydra-genetics create-input-files -d /scratch/wp3/TWIST/${sequencerun} -t N --tc 0 -f && \
hydra-genetics create-input-files -d ${fastqFolder%} -t N --tc 0 -f && \

# Get trio information from SampleSheet.csv into config/samples.tsv
python ${poirotFolder}/extract_sample_sheet_info.py SampleSheet.csv && \

# Start pipeline
snakemake --profile ${poirotFolder}/profiles/slurm/ -s ${poirotFolder}/workflow/Snakefile \
--configfile config/config.yaml -p --latency-wait 5 --rerun-incomplete --keep-going && \

#Copy data back
echo 'Pipeline done and cp back data to OUTBOX' && \
cp /scratch/wp3/TWIST/${sequencerun}/results/multiqc_DNA.html /scratch/wp3/TWIST/${sequencerun}/results/multiqc_${sequencerun}.html
rsync -ru /scratch/wp3/TWIST/${sequencerun}/results/* ${outbox}/${sequencerun}/ && \
rsync -ru /scratch/wp3/TWIST/${sequencerun}/.snakemake ${startDir}/ && \
echo 'All done!'

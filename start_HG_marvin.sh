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
wp3WGS=/projects/wp3/nobackup/WGS

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

rsync -ru SampleSheet.csv /scratch/wp3/TWIST/${sequencerun}/  && \

# Cp data to scratch
echo 'Copy to scratch' && \
rsync -ru SampleSheet.csv /scratch/wp3/TWIST/${sequencerun}/  && \
rsync -ru ${poirotFolder}/config /scratch/wp3/TWIST/${sequencerun}/  && \

cd /scratch/wp3/TWIST/${sequencerun}/  && \

# Create sample and unit files
hydra-genetics create-input-files -d ${fastqFolder} -t N --tc 0 -f && \

# Extract eh WP3 samples from the SampleSheet.csv file
echo 'Create a wp3 only SampleSheet.csv file'   && \
head -n 16 SampleSheet.csv > SampleSheet_WP3.csv   && \
grep WGSWP3 SampleSheet.csv >> SampleSheet_WP3.csv   && \

# extract the WP3 samples from the samples and units files
grep WGSWP3 SampleSheet.csv | cut -f1 -d',' > wp3_samples.txt
head -n 1 samples.tsv  > samples_wp3.tsv  && \
grep -f wp3_samples.txt samples.tsv >> samples_wp3.tsv  && \
head -n 1 units.tsv > units_wp3.tsv  && \
grep -f wp3_samples.txt units.tsv >> units_wp3.tsv   && \
mv samples_wp3.tsv samples.tsv
mv units_wp3.tsv units.tsv

# Add trio and sex info to samples.tsv file and move to config/
python ${poirotFolder}/extract_sample_sheet_info.py SampleSheet_WP3.csv && \

# Start pipeline
snakemake --profile ${poirotFolder}/profiles/slurm/ -s ${poirotFolder}/workflow/Snakefile \
--configfile config/config.yaml -p --latency-wait 5 --rerun-incomplete --keep-going --notemp   && \

#Copy data back
echo 'Pipeline done and cp back data to OUTBOX' && \
cp /scratch/wp3/TWIST/${sequencerun}/results/multiqc_DNA.html /scratch/wp3/TWIST/${sequencerun}/results/multiqc_${sequencerun}.html  && \
rsync -ru /scratch/wp3/TWIST/${sequencerun}/results/* ${outbox}/${sequencerun}/ && \
rsync -ru /scratch/wp3/TWIST/${sequencerun}/.snakemake ${startDir}/ && \
echo 'All done!'

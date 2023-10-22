#!/usr/bin/env bash
eval "$(conda shell.bash hook)"

POIROT_VERSION=bianca_support # version or branch
reference_files=reference_files/ # path to reference files

conda create --name poirot_env python=3.9 -y

conda activate poirot_env

conda install -c conda-forge pip -y

if [ -d poirot_pipeline ];
then
    rm -fr poirot_pipeline
fi

mkdir poirot_pipeline

tar -czvf poirot_pipeline/reference_files.tar.gz ${reference_files}

git clone --branch ${POIROT_VERSION} https://github.com/clinical-genomics-uppsala/poirot_rd_wgs poirot_pipeline/poirot_rd_wgs_${POIROT_VERSION}

pip install -r poirot_pipeline/poirot_rd_wgs_${POIROT_VERSION}/requirements.txt 

# pull the singlarity containers
hydra-genetics prepare-environment create-singularity-files -c config/config.yaml -o poirot_pipeline/singularity_files

conda deactivate

conda-pack -n poirot_env -o poirot_pipeline/poirot_env.tar.gz

mkdir -p poirot_pipeline/hydra-genetics

git clone https://github.com/snakemake/snakemake-wrappers.git poirot_pipeline/snakemake-wrappers

git clone https://github.com/hydra-genetics/alignment.git poirot_pipeline/hydra-genetics/alignment
git clone https://github.com/hydra-genetics/annotation.git poirot_pipeline/hydra-genetics/annotation
git clone https://github.com/hydra-genetics/cnv_sv.git poirot_pipeline/hydra-genetics/cnv_sv
git clone https://github.com/hydra-genetics/compression.git poirot_pipeline/hydra-genetics/compression
git clone https://github.com/hydra-genetics/filtering.git poirot_pipeline/hydra-genetics/filtering
git clone https://github.com/hydra-genetics/misc.git poirot_pipeline/hydra-genetics/misc
git clone https://github.com/hydra-genetics/mitochondrial.git poirot_pipeline/hydra-genetics/mitochondrial
git clone https://github.com/hydra-genetics/prealignment.git poirot_pipeline/hydra-genetics/prealignment
git clone https://github.com/hydra-genetics/parabricks.git poirot_pipeline/hydra-genetics/parabricks
git clone https://github.com/hydra-genetics/qc.git poirot_pipeline/hydra-genetics/qc
git clone https://github.com/hydra-genetics/snv_indels.git poirot_pipeline/hydra-genetics/snv_indels


tar -zcvf poirot_pipeline.tar.gz poirot_pipeline

if [ -d poirot_pipeline ];
then
    rm -fr poirot_pipeline
fi
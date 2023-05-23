#!/usr/bin/env bash
eval "$(conda shell.bash hook)"

POIROT_VERSION=develop

conda create --name poirot_rd_wgs_${POIROT_VERSION} python=3.9 -y

conda activate poirot_rd_wgs_${POIROT_VERSION}

conda install -c conda-forge pip -y

if [ -d poirot_rd_wgs_${POIROT_VERSION} ];
then
    rm -fr poirot_rd_wgs_${POIROT_VERSION}
fi

mkdir poirot_rd_wgs_${POIROT_VERSION}

git clone --branch ${POIROT_VERSION} https://github.com/clinical-genomics-uppsala/poirot_rd_wgs poirot_rd_wgs_${POIROT_VERSION}/poirot_rd_wgs

pip install -r poirot_rd_wgs_${POIROT_VERSION}/poirot_rd_wgs/requirements.txt 

conda deactivate

conda-pack -n poirot_rd_wgs_${POIROT_VERSION} -o poirot_rd_wgs_${POIROT_VERSION}/env.tar.gz

mkdir -p poirot_rd_wgs_${POIROT_VERSION}/hydra-genetics

git clone https://github.com/snakemake/snakemake-wrappers.git poirot_rd_wgs_${POIROT_VERSION}/snakemake-wrappers


git clone https://github.com/hydra-genetics/alignment.git poirot_rd_wgs_${POIROT_VERSION}/hydra-genetics/alignment
git clone https://github.com/hydra-genetics/annotation.git poirot_rd_wgs_${POIROT_VERSION}/hydra-genetics/annotation
git clone https://github.com/hydra-genetics/cnv_sv.git poirot_rd_wgs_${POIROT_VERSION}/hydra-genetics/cnv_sv
git clone https://github.com/hydra-genetics/compression.git poirot_rd_wgs_${POIROT_VERSION}/hydra-genetics/compression
git clone https://github.com/hydra-genetics/filtering.git poirot_rd_wgs_${POIROT_VERSION}/hydra-genetics/filtering
git clone https://github.com/hydra-genetics/misc.git poirot_rd_wgs_${POIROT_VERSION}/hydra-genetics/misc
git clone https://github.com/hydra-genetics/mitochondrial.git poirot_rd_wgs_${POIROT_VERSION}/hydra-genetics/mitochondrial
git clone https://github.com/hydra-genetics/prealignment.git poirot_rd_wgs_${POIROT_VERSION}/hydra-genetics/prealignment
git clone https://github.com/hydra-genetics/parabricks.git poirot_rd_wgs_${POIROT_VERSION}/hydra-genetics/parabricks
git clone https://github.com/hydra-genetics/qc.git poirot_rd_wgs_${POIROT_VERSION}/hydra-genetics/qc
git clone https://github.com/hydra-genetics/snv_indels.git poirot_rd_wgs_${POIROT_VERSION}/hydra-genetics/snv_indels


tar -zcvf poirot_rd_wgs_${POIROT_VERSION}.tar.gz poirot_rd_wgs_${POIROT_VERSION}

if [ -d poirot_rd_wgs_${POIROT_VERSION} ];
then
    rm -fr poirot_rd_wgs_${POIROT_VERSION}
fi
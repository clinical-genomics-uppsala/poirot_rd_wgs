# Packaging Poirot for offline environments

The build_conda.sh script packages the pipeline, cofigs, reference files and apptainer files for an offline environment. 

When succesfully run the script will generate four compressed tar archives:

* design_and_ref_files.tar.gz
* poirot_config_${PIPELINE_VERSION}.tar.gz
* apptainer_cache.tar.gz
* poirot_config_${CONFIG_VERSION}.tar.gz

The requirments listed in requirements.txt are packaged using conda-pack in a .tar.gz in the poirot_config_${PIPELINE_VERSION}.tar.gz. The snakemake-wrappers github repo and all hydra-genetics modules required by the pipeline are cloned and packaged in poirot_config_${PIPELINE_VERSION}.tar.gz.

```
export TAG_OR_BRANCH="v0.11.0"
export CONFIG_VERSION="v0.14.0"
export PIPELINE_NAME="poirot_rd_wgs"
export PYTHON_VERSION="3.9"
export PIPELINE_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/poirot_rd_wgs.git"
export CONFIG_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/poirot_config.git"

```

Download only pipeline files
```bash
bash build_conda.sh --pipeline-only 
```

Download only config files
```bash
bash build_conda.sh --config-only
```

Download only singularity files
```bash
bash build_conda.sh --containers-only
```

Download only reference and design files

```bash
bash build_conda.sh --references-only poirot_config/config/references/gene_panels.hg38.yaml  poirot_config/config/references/mito_refs.hg38.yaml poirot_config/config/references/references.hg38.yaml poirot_config/config/references/str_panels.hg38.yaml melt_refs.hg38.yaml vep.hg38.yaml
```

Download only gene panels

```bash
bash build_conda.sh --references-only poirot_config/config/references/gene_panels.hg38.yaml  
```

Download all files

```bash
bash build_conda.sh --all poirot_config/config/references/gene_panels.hg38.yaml  poirot_config/config/references/mito_refs.hg38.yaml poirot_config/config/references/references.hg38.yaml poirot_config/config/references/str_panels.hg38.yaml melt_refs.hg38.yaml vep.hg38.yaml
```
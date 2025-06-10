#Building version v0.8.0

```bash
TAG_OR_BRANCH="v0.8.0" CONFIG_VERSION="v0.10.0" PIPELINE_NAME="poirot_rd_wgs" PYTHON_VERSION="3.9" \
PIPELINE_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/poirot_rd_wgs.git" \
CONFIG_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/poirot_config.git" \
bash build_conda.sh poirot_config/config/references/gene_panels.hg38.yaml  poirot_config/config/references/mito_refs.hg38.yaml \
poirot_config/config/references/references.hg38.yaml poirot_config/config/references/str_panels.hg38.yaml
```
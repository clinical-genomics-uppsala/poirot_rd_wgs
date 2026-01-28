#!/usr/bin/env bash
set -e

eval "$(conda shell.bash hook)"

# Function to download and package the pipeline
download_pipeline() {
    echo "=== Downloading Pipeline ==="
    
    # Create and activate conda environment in the current directory, then install pipeline requirements
    echo "Creating conda environment: ${PIPELINE_NAME}_${TAG_OR_BRANCH}_env"
    mamba create --prefix ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env python=${PYTHON_VERSION} -y
    conda activate ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env
    
    # Clean up existing pipeline directory if it exists
    if [ -d ${PIPELINE_NAME}_${TAG_OR_BRANCH} ]; then
        echo "Removing existing pipeline directory: ${PIPELINE_NAME}_${TAG_OR_BRANCH}"
        rm -fr ${PIPELINE_NAME}_${TAG_OR_BRANCH}
    fi
    
    mkdir ${PIPELINE_NAME}_${TAG_OR_BRANCH}
    
    # Clone the required version of the pipeline
    echo "Cloning pipeline from ${PIPELINE_GITHUB_REPO} (branch: ${TAG_OR_BRANCH})"
    git clone --branch ${TAG_OR_BRANCH} ${PIPELINE_GITHUB_REPO} ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}
    
    # Install the requirements for the pipeline
    echo "Installing pipeline requirements"
    export PYTHONNOUSERSITE=1 # stops pip looking in ˙$HOME/.local for packages
    ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env/bin/pip3 install -r ${PIPELINE_NAME}_${TAG_OR_BRANCH}/${PIPELINE_NAME}/requirements.txt
    
    # Pack the environment with the requirements installed
    echo "Packing conda environment"
    conda pack --prefix ./${PIPELINE_NAME}_${TAG_OR_BRANCH}_env -o ${PIPELINE_NAME}_${TAG_OR_BRANCH}/env.tar.gz
    
    # Clone snakemake-wrappers and hydra-genetics modules
    echo "Cloning snakemake-wrappers and hydra-genetics modules"
    mkdir -p ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics
    
    git clone https://github.com/snakemake/snakemake-wrappers.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/snakemake-wrappers
    
    # Array of hydra-genetics modules to clone
    
    local hydra_modules=(
        "alignment"
        "annotation" 
        "biomarker"
        "cnv_sv"
        "compression"
        "filtering"
        "fusions"
        "misc"
        "mitochondrial"
        "parabricks"
        "prealignment"
        "qc"
        "reports"
        "snv_indels"
        "references"
    )
    
    # Clone each hydra-genetics module
    for module in "${hydra_modules[@]}"; do
        echo "Cloning hydra-genetics/${module}"
        git clone https://github.com/hydra-genetics/${module}.git ${PIPELINE_NAME}_${TAG_OR_BRANCH}/hydra-genetics/${module}
    done
    
    # Pack all cloned repositories
    echo "Creating pipeline archive: ${PIPELINE_NAME}_${TAG_OR_BRANCH}.tar.gz"
    tar -zcvf ${PIPELINE_NAME}_${TAG_OR_BRANCH}.tar.gz ${PIPELINE_NAME}_${TAG_OR_BRANCH}
    
    # Download the config files from the config repo
    echo "Downloading config files from ${CONFIG_GITHUB_REPO} (version: ${CONFIG_VERSION})"
    git clone --branch ${CONFIG_VERSION} ${CONFIG_GITHUB_REPO} poirot_config/
    
    echo "Pipeline download completed successfully"
}

# Function to download containers
download_containers() {
    echo "=== Downloading Containers ==="
    
    # Check if config directory exists, if not download it
    if [ ! -d poirot_config ]; then
        echo "Config directory not found, downloading config files"
        git clone --branch ${CONFIG_VERSION} ${CONFIG_GITHUB_REPO} poirot_config/
    fi
    
    # Download containers using hydra-genetics
    echo "Creating singularity files using hydra-genetics"
    hydra-genetics prepare-environment create-singularity-files -c poirot_config/config/config_production_pipeline.yaml -o apptainer_cache
    
    # Copy additional container (MELT)
    echo "Copying MELT container"
    cp /projects/wp3/Software/MELTv2.2.2/MELT_v2.2.2.sif apptainer_cache
    
    # Create container archive
    echo "Creating container archive: apptainer_cache.tar.gz"
    tar -czvf apptainer_cache.tar.gz apptainer_cache
    
    echo "Container download completed successfully"
}

# Function to download design and reference files
download_design_and_reference_files() {
    echo "=== Downloading Design and Reference Files ==="
    
    # Check if we have any reference configs to process
    if [ $# -eq 0 ]; then
        echo "No reference config files provided"
        return
    fi
    
    # Check if config directory exists, if not download it
    if [ ! -d poirot_config ]; then
        echo "Config directory not found, downloading config files"
        git clone --branch ${CONFIG_VERSION} ${CONFIG_GITHUB_REPO} poirot_config/
    fi
    
    # Download references for each provided config file
    for reference_config in "$@"; do
        echo "Processing reference config: ${reference_config}"
        if [ -f "$reference_config" ]; then
            hydra-genetics --debug references download -o design_and_ref_files -v "$reference_config"
        else
            echo "Warning: Reference config file not found: ${reference_config}"
        fi
    done
    
    # Create reference files archive only if directory exists and is not empty
    if [ -d design_and_ref_files ] && [ "$(ls -A design_and_ref_files 2>/dev/null)" ]; then
        echo "Creating reference files archive: design_and_ref_files.tar.gz"
        tar -czvf design_and_ref_files.tar.gz design_and_ref_files
        echo "Design and reference files download completed successfully"
    else
        echo "No reference files were downloaded"
    fi
}

# Function to download only config files
download_config() {
    echo "=== Downloading Config Files ==="
    
    # Clean up existing config directory if it exists
    if [ -d poirot_config ]; then
        echo "Removing existing config directory: poirot_config"
        rm -fr poirot_config
    fi
    
    # Download the config files from the config repo
    echo "Downloading config files from ${CONFIG_GITHUB_REPO} (version: ${CONFIG_VERSION})"
    git clone --branch ${CONFIG_VERSION} ${CONFIG_GITHUB_REPO} poirot_config_${CONFIG_VERSION}
    
    ## add the pipeline version to the profiles config files
    envsubst < poirot_config_${CONFIG_VERSION}/profiles/${PROFILE_NAME}/config.yaml > poirot_config_${CONFIG_VERSION}/profiles/${PROFILE_NAME}/config.yaml.sub
    mv poirot_config_${CONFIG_VERSION}/profiles/${PROFILE_NAME}/config.yaml.sub poirot_config_${CONFIG_VERSION}/profiles/${PROFILE_NAME}/config.yaml

    # Create config archive with version number
    echo "Creating config archive: poirot_config_${CONFIG_VERSION}.tar.gz"
    tar -czvf poirot_config_${CONFIG_VERSION}.tar.gz poirot_config_${CONFIG_VERSION}
    
    echo "Config files download completed successfully"
}

# Function to clean up temporary directories and files
cleanup() {
    echo "=== Cleaning up temporary files ==="
    
    # Clean pipeline-related files
    if [ "$DOWNLOAD_PIPELINE" = true ]; then
        if [ -d ${PIPELINE_NAME}_${TAG_OR_BRANCH}_env ]; then
            echo "Removing temporary conda environment: ${PIPELINE_NAME}_${TAG_OR_BRANCH}_env"
            rm -fr ${PIPELINE_NAME}_${TAG_OR_BRANCH}_env
        fi
        
        if [ -d ${PIPELINE_NAME}_${TAG_OR_BRANCH} ]; then
            echo "Removing temporary pipeline directory: ${PIPELINE_NAME}_${TAG_OR_BRANCH}"
            rm -fr ${PIPELINE_NAME}_${TAG_OR_BRANCH}
        fi
    fi
    
    # Clean config directory (shared between containers, references, and config-only)
    if [ -d poirot_config ]; then
        echo "Removing temporary config directory: poirot_config"
        rm -fr poirot_config
    fi
    
    # Clean container-related files
    if [ "$DOWNLOAD_CONTAINERS" = true ]; then
        if [ -d apptainer_cache ]; then
            echo "Removing temporary container cache: apptainer_cache"
            rm -fr apptainer_cache
        fi
    fi
    
    # Clean reference-related files  
    if [ "$DOWNLOAD_REFERENCES" = true ]; then
        if [ -d design_and_ref_files ]; then
            echo "Removing temporary reference files directory: design_and_ref_files"
            rm -fr design_and_ref_files
        fi
    fi
    
    echo "Cleanup completed"
}

# Function to validate required environment variables
validate_environment() {
    local required_vars=("TAG_OR_BRANCH" "CONFIG_VERSION" "PIPELINE_NAME" "PYTHON_VERSION" "PIPELINE_GITHUB_REPO" "CONFIG_GITHUB_REPO")
    local missing_vars=()
    
    for var in "${required_vars[@]}"; do
        if [ -z "${!var}" ]; then
            missing_vars+=("$var")
        fi
    done
    
    if [ ${#missing_vars[@]} -gt 0 ]; then
        echo "Error: Missing required environment variables:"
        printf "  - %s\n" "${missing_vars[@]}"
        echo ""
        echo "Please set all required variables before running this script."
        echo "Example usage:"
        echo 'TAG_OR_BRANCH="v0.8.0" CONFIG_VERSION="v0.10.0" PIPELINE_NAME="poirot_rd_wgs" \\'
        echo 'PYTHON_VERSION="3.9" PIPELINE_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/poirot_rd_wgs.git" \\'
        echo 'CONFIG_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/poirot_config.git" \\'
        echo 'bash build_conda.sh config1.yaml config2.yaml ...'
        exit 1
    fi
}

# Function to display usage information
usage() {
    cat << EOF
Usage: $0 [OPTIONS] [reference_config1.yaml] [reference_config2.yaml] ...

This script builds a complete pipeline package including:
  1. Pipeline code and conda environment
  2. Container images (Singularity/Apptainer)
  3. Design and reference files

OPTIONS:
  -h, --help              Show this help message and exit
  -p, --pipeline-only     Download only the pipeline (step 1)
  -c, --containers-only   Download only the containers (step 2)
  -r, --references-only   Download only the design and reference files (step 3)
  -g, --config-only       Download only the config files
  -a, --all              Download all components (default behavior)

If no options are specified, all components will be downloaded.

Required environment variables:
  - TAG_OR_BRANCH: Git tag or branch of the pipeline to build
  - CONFIG_VERSION: Version of the configuration repository
  - PIPELINE_NAME: Name of the pipeline
  - PYTHON_VERSION: Python version for conda environment
  - PIPELINE_GITHUB_REPO: URL of the pipeline repository
  - CONFIG_GITHUB_REPO: URL of the configuration repository

Examples:
  # Download all components (default)
  TAG_OR_BRANCH="v0.8.0" CONFIG_VERSION="v0.10.0" PIPELINE_NAME="poirot_rd_wgs" \\
  PYTHON_VERSION="3.9" \\
  PIPELINE_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/poirot_rd_wgs.git" \\
  CONFIG_GITHUB_REPO="https://github.com/clinical-genomics-uppsala/poirot_config.git" \\
  bash $0 poirot_config/config/references/gene_panels.hg38.yaml

  # Download only the pipeline
  bash $0 --pipeline-only

  # Download only containers
  bash $0 --containers-only

  # Download only config files
  bash $0 --config-only

  # Download only reference files
  bash $0 --references-only poirot_config/config/references/gene_panels.hg38.yaml \\
                            poirot_config/config/references/mito_refs.hg38.yaml

EOF
}

# Function to parse command line arguments
parse_arguments() {
    DOWNLOAD_PIPELINE=false
    DOWNLOAD_CONTAINERS=false
    DOWNLOAD_REFERENCES=false
    DOWNLOAD_CONFIG=false
    REFERENCE_CONFIGS=()
    
    # If no arguments, download all components (except config-only)
    if [ $# -eq 0 ]; then
        DOWNLOAD_PIPELINE=true
        DOWNLOAD_CONTAINERS=true
        DOWNLOAD_REFERENCES=true
        return
    fi
    
    while [[ $# -gt 0 ]]; do
        case $1 in
            -h|--help)
                usage
                exit 0
                ;;
            -p|--pipeline-only)
                DOWNLOAD_PIPELINE=true
                shift
                ;;
            -c|--containers-only)
                DOWNLOAD_CONTAINERS=true
                shift
                ;;
            -r|--references-only)
                DOWNLOAD_REFERENCES=true
                shift
                ;;
            -g|--config-only)
                DOWNLOAD_CONFIG=true
                shift
                ;;
            -a|--all)
                DOWNLOAD_PIPELINE=true
                DOWNLOAD_CONTAINERS=true
                DOWNLOAD_REFERENCES=true
                shift
                ;;
            *.yaml|*.yml)
                REFERENCE_CONFIGS+=("$1")
                shift
                ;;
            *)
                echo "Unknown option: $1"
                usage
                exit 1
                ;;
        esac
    done
    
    # If no specific component selected, download all (except config-only)
    if [ "$DOWNLOAD_PIPELINE" = false ] && [ "$DOWNLOAD_CONTAINERS" = false ] && [ "$DOWNLOAD_REFERENCES" = false ] && [ "$DOWNLOAD_CONFIG" = false ]; then
        DOWNLOAD_PIPELINE=true
        DOWNLOAD_CONTAINERS=true
        DOWNLOAD_REFERENCES=true
    fi
    
    # If downloading references but no config files provided, warn user
    if [ "$DOWNLOAD_REFERENCES" = true ] && [ ${#REFERENCE_CONFIGS[@]} -eq 0 ]; then
        echo "Warning: References download requested but no reference config files provided"
        echo "Reference files will not be downloaded without config files"
        DOWNLOAD_REFERENCES=false
    fi
}

# Main execution function
main() {
    # Parse command line arguments
    parse_arguments "$@"
    
    # Validate environment variables
    validate_environment
    
    echo "Starting build process for ${PIPELINE_NAME} ${TAG_OR_BRANCH}"
    echo "Python version: ${PYTHON_VERSION}"
    echo "Config version: ${CONFIG_VERSION}"
    echo "Pipeline repo: ${PIPELINE_GITHUB_REPO}"
    echo "Config repo: ${CONFIG_GITHUB_REPO}"
    echo ""
    echo "Components to download:"
    echo "  - Pipeline: $DOWNLOAD_PIPELINE"
    echo "  - Containers: $DOWNLOAD_CONTAINERS"  
    echo "  - References: $DOWNLOAD_REFERENCES"
    echo "  - Config: $DOWNLOAD_CONFIG"
    if [ ${#REFERENCE_CONFIGS[@]} -gt 0 ]; then
        echo "  - Reference configs: ${REFERENCE_CONFIGS[*]}"
    fi
    echo ""
    
    # Execute selected build steps
    if [ "$DOWNLOAD_PIPELINE" = true ]; then
        download_pipeline
    fi
    
    if [ "$DOWNLOAD_CONTAINERS" = true ]; then
        download_containers
    fi
    
    if [ "$DOWNLOAD_REFERENCES" = true ]; then
        download_design_and_reference_files "${REFERENCE_CONFIGS[@]}"
    fi
    
    if [ "$DOWNLOAD_CONFIG" = true ]; then
        download_config
    fi
    
    # Clean up temporary files
    cleanup
    
    # Deactivate conda environment if it was activated
    if [ "$DOWNLOAD_PIPELINE" = true ]; then
        conda deactivate
    fi
    
    echo ""
    echo "Build process completed successfully!"
    echo "Generated files:"
    if [ "$DOWNLOAD_PIPELINE" = true ]; then
        echo "  - ${PIPELINE_NAME}_${TAG_OR_BRANCH}.tar.gz (pipeline)"
    fi
    if [ "$DOWNLOAD_CONTAINERS" = true ]; then
        echo "  - apptainer_cache.tar.gz (containers)"
    fi
    if [ "$DOWNLOAD_REFERENCES" = true ] && [ ${#REFERENCE_CONFIGS[@]} -gt 0 ]; then
        echo "  - design_and_ref_files.tar.gz (reference files)"
    fi
    if [ "$DOWNLOAD_CONFIG" = true ]; then
        echo "  - poirot_config_${CONFIG_VERSION}.tar.gz (config files)"
    fi
}

# Run main function with all provided arguments
main "$@"
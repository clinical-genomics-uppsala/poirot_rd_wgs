# Running the pipeline

## Requirements
**Recommended hardware**

- CPU: >10 cores per sample
- Memory: 6GB per core
- Storage: >75GB per sample

**Note**: Running the pipeline with less resources may work, but has not been tested.

**Software**

- [python](https://www.python.org/), version 3.9 or newer
- [pip3](https://pypi.org/project/pip/)
- [virtuelenv](https://docs.python.org/3/library/venv.html)
- [singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html)

**Nice to have**

- DRMAA compatible scheduler


## Installation
A list of releases of Poirot can be found at: [Releases](https://github.com/clinical-genomics-uppsala/poirot_rd_wgs/releases).

### Clone the Poirot git repo
We recommend that the repository is cloned to your working directory. 
```bash
# Set up a working directory path
WORKING_DIRECTORY="/path_working_to_directory"
```

Fetch pipeline
```bash
# Set version
VERSION="v0.5.1"

# Clone selected version
git clone --branch ${VERSION} https://github.com/clinical-genomics-uppsala/poirot_rd_wgs.git ${WORKING_DIRECTORY}
```

### Create python environment
To run the Poirot pipeline a python virtual environment is needed.
```bash
# Enter working directory
cd ${WORKING_DIRECTORY}

# Create a new virtual environment
python3 -m venv ${WORKING_DIRECTORY}/virtual/environment
```

### Install pipeline requirements
Activate the virtual environment and install pipeline requirements specified in `requirements.txt`.
```bash
# Enter working directory
cd ${WORKING_DIRECTORY}

# Activate python environment
source environment/bin/activate

# Install requirements
pip install -r requirements.txt
```

## Input sample files
The pipeline uses sample input files (`samples.tsv` and `units.tsv`) with information regarding sample information, sequencing meta information as well as the location of the fastq-files. Specification for the input files can be found at [Poirot schemas](https://github.com/clinical-genomics-uppsala/poirot_rd_wgs/tree/main/workflow/schemas). Using the python virtual environment created above it is possible to generate these files automatically using [hydra-genetics create-input-files](https://hydra-genetics.readthedocs.io/en/latest/run_pipeline/create_sample_files/):
```bash
hydra-genetics create-input-files -d path/to/fastq-files/
```

## Run command
Using the activated python virtual environment created above, this is a basic command for running the pipeline:
```bash
snakemake --profile profiles/NAME_OF_PROFILE -s workflow/Snakefile
```  
<br />
The are many additional [snakemake running options](https://snakemake.readthedocs.io/en/stable/executing/cli.html#) some of which is listed below. However, options that are always used should be put in the [profile](https://hydra-genetics.readthedocs.io/en/latest/run_pipeline/profile/).

* --notemp - Saves all intermediate files. Good for development and testing different options.
* --until <rule> - Runs only rules dependent on the specified rule.

<br />
**Note:** Remember to have singularity and drmaa available on the system where the pipeline will be run.

<br />
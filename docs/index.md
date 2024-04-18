# Welcome to Poirot

<p align="center">
<a href="https://github.com/clinical-genomics-uppsala/poirot_rd_wgs/">https://github.com/clinical-genomics-uppsala/poirot_rd_wgs/</a>
</p>

This pipeline is run at Clinical Genomics Uppsala to call variants from short-read illumina WGS data from rare disease patients. This [snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline is built using module system from [Hydra Genetics](https://github.com/hydra-genetics/) to process paired-end `.fastq.gz` files from Illumina whole genome sequencing.

Poirot uses the following hydra genetics modules:

- [Alignment](https://github.com/hydra-genetics/alignment/tree/1c54479)
- [annotation](https://github.com/hydra-genetics/alignment/tree/v0.3.0)
- [compression](https://github.com/hydra-genetics/alignment/tree/v2.0.0)
- [cnv_sv](https://github.com/hydra-genetics/alignment/tree/v0.5.0)
- [filtering](https://github.com/hydra-genetics/alignment/tree/v0.3.0)
- [parabricks](https://github.com/hydra-genetics/alignment/tree/1.2.0)
- [prealignment](https://github.com/hydra-genetics/alignment/tree/v1.2.0)
- [qc](https://github.com/hydra-genetics/alignment/tree/da66130)
- [mitochondrial](https://github.com/hydra-genetics/alignment/tree/v0.1.0)
- [snv_indels](https://github.com/hydra-genetics/alignment/tree/3935ecf)

# Hydra-genetics

We are an organization/community with the goal of making [snakemake](https://snakemake.readthedocs.io/en/stable/index.html) pipeline development easier, faster, a bit more structured and of higher quality.

We do this by providing [snakemake modules](https://snakemake.readthedocs.io/en/stable/snakefiles/modularization.html#modules) that can be combined to create a complete analysis or included in already existing pipelines. All modules are subjected to extensive testing to make sure that new releases doesn't unexpectedly break existing pipeline or deviate from guidelines and best practices on how to write code.

# Snakemake
Poirot and Hydra-genetics are snakemake bases pipeline/tools. The [Snakemake](https://snakemake.readthedocs.io/en/stable/index.html) workflow management system is a tool to create reproducible and scalable data analyses. Workflows are described via a human readable, Python based language. They can be seamlessly scaled to server, cluster, grid and cloud environments, without the need to modify the workflow definition. Finally, Snakemake workflows can entail a description of required software, which will be automatically deployed to any execution environment. 

If Snakemake is new to you a good place to start is doing the [snakemake tutorial](https://snakemake.readthedocs.io/en/stable/tutorial/tutorial.html) since this will help you setting Poirot up.
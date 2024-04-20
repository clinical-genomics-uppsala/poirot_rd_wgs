# Changelog

### [0.4.1](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/compare/v0.4.0...v0.4.1) (2024-04-20)


### Bug Fixes

* handle col creation when no trios present ([6f55abf](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/6f55abf57cef25d8c438e6bdf2d2c621498d6b6a))
* **Snakefile:** add rule order for bgzip and rename rule ([60d4ee8](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/60d4ee837e2e17669c1c91526365f034916d4e40))


### Documentation

* Update .readthedocs.yaml ([c740d3b](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/c740d3b15184b346d04467170a74268d90b67d9e))
* Update README.md ([871dc9d](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/871dc9d620ca2098a7b96dcb7495307a8bc712ec))

## [0.4.0](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/compare/v0.3.1...v0.4.0) (2024-04-18)


### Features

* add peddy html report to output_list.json ([949e5a2](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/949e5a2b27a1e79f017149c2e77915592498954b))
* improve handling of config validation ([603b258](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/603b258aa3ad486c982c3196aa2992500243108e))
* remove temp on peddy sex check file ([7b16868](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/7b168685a759adad9058e296254864fc27ac4a71))
* update hydra-modules ([44857c8](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/44857c82031dff122cad72258574471751bcd385))
* update the alignment module version ([9e58777](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/9e587770af185c91af67a99f91beab36a2ddd9f2))
* use misc module for bgzip and tabix ([07f21c2](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/07f21c21eb60041638dcf73027248c8d20e181e5))
* use run_deepvariant script in updated snv_indels module ([6d4dfc2](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/6d4dfc2a30ae997c99b9dd695a0102ed9cf7c3e9))


### Bug Fixes

* update manta input names ([843e4f0](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/843e4f0c22bc53415ad53bec27296c3baf4306b4))


### Documentation

* start using rtd ([9a6bc31](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/9a6bc31fd08981c8820e904f2d6d835d857824de))
* update dag graph ([42d6c12](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/42d6c12858c3d0d8bba8feba0393b3e628918739))
* update dag image ([2aa37f1](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/2aa37f1f9be78024a0268e69574810d8a14040b9))
* update docs ([ddf7a5e](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/ddf7a5e451bec470985931c984c84216c050ca34))
* update mkdocs nav ([4df2996](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/4df29962e5092f6f94411d013024575f9b535a29))
* update README ([f9544c3](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/f9544c32b8d19a85be78b792055fcdd876d86341))
* update REAME ([3737764](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/3737764818ce9a2851e691d4f914f8266f6d87e7))
* update rulegraph ([2e6f88b](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/2e6f88bf9da9e587691cbdde5c1226b770c41653))

### [0.3.1](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/compare/v0.3.0...v0.3.1) (2024-02-28)


### Bug Fixes

* handle incomplete trios ([35e7e5a](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/35e7e5a48e27932d3f13f60e2b72ca1da988f87a))

## [0.3.0](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/compare/v0.2.0...v0.3.0) (2024-01-03)


### Features

* add verifybamid2 to the QC ([540b6d6](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/540b6d66d996f0cd6bd3947bfc6c9335e19d09ba))
* handle when the coverage bed file that has mitochondria genes ([75d6296](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/75d62969a51f29debf059933d7259f001d6bfa90))

## [0.2.0](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/compare/v0.1.0...v0.2.0) (2023-12-07)


### Features

* add constraint to the rule resources ([bc7fd42](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/bc7fd42aa83f86cceed0f7cbd957e9a1828e0497))
* add singularity download and reference files compression ([af55a2a](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/af55a2aabc80d9422fa01447868f9a906a224611))
* add snakemake wrapper prefix path to bianca profile ([ed90e49](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/ed90e490be35b853ab0b05f328c761480324fb96))
* change deeptrio to deepvariant in the upd analysis ([10d93a7](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/10d93a7044338e08c59a5bf0a08c51c89d318a6e))
* support running the pipeline on bianca cluster ([1691c4d](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/1691c4d0b12ec2e0bc3c565b0c44950f32ca6c5a))


### Bug Fixes

* point to the WGS model for deepvariant CPU ([18c4966](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/18c4966edba7838c899ed7e4710e992b1804b91c))
* use updated qc module so that samtools stats avoids using design_bed ([c970216](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/c9702168878feaedfaf3e4a1de536380b3505d4a))

## 0.1.0 (2023-10-20)


### Features

* add conversion from stranger vcf to bed file ([c0e8c79](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/c0e8c79661f26d8c4a410b02a85ea2ab5d44d99e))
* add deeptrio and upd ([38ddc8f](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/38ddc8f40202262f206802942ce46570c63032b3))
* Add ExpansionHunter ([43b0c50](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/43b0c509470f0b3489a7ff19e1ddf6dcf12fe93b))
* add functions to generate cp rules ([fdd4ce9](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/fdd4ce9c84da7f35b6ff3ec467867a6291ea4472))
* add logging to coverage script ([80d5244](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/80d5244ad2c27783f670094068e03def6aceefff))
* Add options to run the bwa mapping  and deepvariant calling on CPU ([6493c85](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/6493c850ac175baeffc4f050a66c971e01d833e1))
* add reference to VCF header of tiddit and svdb merge VCF files ([d301d04](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/d301d04d62ddfca29f58e620f438e0bffa3e8a9e))
* add scripts for sample info extraction and running the pipeline ([8a0530c](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/8a0530c563c4b2255b78d755bb655071bc3e8ad8))
* add SMNcopynumbercaller ([7fcdc94](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/7fcdc94997600a1a41347b26fad63bf509f03c1d))
* add stranger rule ([af4e15a](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/af4e15ad7307c0a8be6891dc8b58d5ee3fa18805))
* add vcf to aed for cnvpytor vcf files ([75e7970](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/75e79703e958a37efafcaba79252166ea84a1cc2))
* change how the reference is added in the header of  the vcfs ([32ab5b1](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/32ab5b12605c772ec826094fe80370c429a8cb93))
* change to samtools cram ([65dd919](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/65dd9192acd5e2ddc507cbac03ce32432d84c0b6))
* handle cases in the sample sheet where the sex is unknown ([ad0c5b1](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/ad0c5b1478b17a44748b73d3f432d683736596a9))
* handle update gres requirements ([8807e2c](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/8807e2c407b9000a3c124b012b454b56c0aff86b))
* increase the runtime for cnvpytor ([e571362](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/e57136210ed21e9f6b4d557786d0e86ee0d2405c))
* look for sample order and replacement files in config folder ([a0b67de](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/a0b67ded9a33ae863b99f9b82317ed655cc1410b))
* make deepvariant vcf compressed to work with script ([d00ab67](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/d00ab67bf8e0d3abe19f41e7eb9b8be28bb49b92))
* order samples in multiqc based on sample sheet ([a8d5ac1](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/a8d5ac1c3505b71b7fd917384bb949f80ebe8f83))
* remove rule and script for sample ordering ([5117352](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/5117352d6d308e759ae65ddcd82010113b37a777))
* update cnvpytor version ([4574806](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/4574806e2bbb3fd95522b3311ae7c786a8747156))
* update create_peddy_fam.py to use sample.tsv ([c48f52d](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/c48f52d9b11739b7ddfaa9cf62d4d8236bfafbd8))
* update GPU partition info ([11e7eb8](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/11e7eb8c724e0b841933524b6fd05fc219a0a07e))
* update marvin start script to process only WP3 samples in the sample sheet ([7d1ea28](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/7d1ea28f2368e79120ee13fb22a707cc0a1d0622))
* update parabricks steps for parabricks 4 ([02eda9b](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/02eda9bc246e7056dee9dbe699b3c7b04a894370))
* update reference file ([622f296](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/622f296cdf6b7f86cc7afc8c1944b50862105c0c))
* Update start_HG_marvin.sh ([67571f6](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/67571f61b9f191d2b0fae4a205d71093109ad1d2))
* update the config files ([3003d01](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/3003d017f3e2e2583a7c1765cb5fff8c7c6eff6b))
* update to latest versioned hydra-genetics modules ([3074822](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/30748225a0bcce72898d0fa1bed1d1355c94b71e))


### Bug Fixes

* **_cp_reviewer:** set directory on output of rule ([3b69884](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/3b698849629dc52938897f28fd875c2580f3a153))
* add /projects,/scratch,/data to --bind in profile ([16e90df](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/16e90df8d82d695a24ac4ff0cdf7fe57226e20f4))
* add bam index a input for CNVpytor ([444bbb3](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/444bbb36cce7a254cc26c9b68069694b1abff652))
* add bam index for child in deeptrio make_examples ([7303275](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/730327549c300b4878097fa719f2c1cb5adc7446))
* add configfile to multiqc ([ef7d36a](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/ef7d36afe473a0bf154731fbc759456669a2fb3b))
* Change how gene panel files paths are constructed ([e0a84f7](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/e0a84f7ca7a6611fa3d4335cdc33835499dbf2d5))
* change to uncompressed bed file ([246415e](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/246415e44a7d23ff9a011a617a846c72188202b1))
* change wildcard constraints to get the _copy_spring rule working ([df216c9](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/df216c95fa213d2a972f532481810304b323a4bf))
* correctly name family id column ([c2c838e](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/c2c838ee57e6b546cf87f9b68ec767063a90e824))
* directory and vep vonfig ([2ec1563](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/2ec156343d4150a70fa6863930441a23176ad58b))
* fix multiqc params.extra ([0ab3279](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/0ab3279f8c7d4adc189920fe9647ee2ab331209b))
* fix path for log and benchmark rules ([93e60d6](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/93e60d6d74f82ac342c38f2ad5b61a9f52cdb75c))
* fix pbrun deepvariant command for v1.0.0 module ([c76146a](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/c76146a14c94e10617ceaeb854911c93abf88c33))
* handle cases when all members of a trio are not in the samples.tsv ([1394fdc](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/1394fdc8cdea88271ca44af85a45e5f40e330fb7))
* Handle cases where the rel_df is empty due to lack of trios or peddy errors ([08187a8](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/08187a8abc10cb40d674bf3e224b6d81b587d040))
* handle overlapping genes by sorting by gene name ([dd7497d](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/dd7497dade15334263203f1b763fb2c0014f4ff6))
* handle overlapping genes in create_excel.py ([91f09e7](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/91f09e7cc1d807656f8e32a665723e14305712e8))
* Include all peddy files as inputs for multiqc ([76113ab](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/76113ab082a74b840268f4a9600fb5d37f3e553c))
* keep the gVCF, useful when runs crash ([955dc8d](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/955dc8db4eacb1bf4995a7e34ba63d1496a625d1))
* list all input files in create_cov_excel rile ([53171aa](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/53171aa7512d85eccde14d14614165b5d3039d36))
* make output_file for _copy_reviewer a directory ([07f1e94](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/07f1e94f35b3d5bbb6a41da42c2fcd4f6508501e))
* **peddy:** point to the samples in config for input ([bd217f6](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/bd217f615a691be0e868aae6a4766fa27df2f15e))
* point to the correct module versions ([5b124ba](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/5b124bab06b365a6b6888955a9de565535b0ee61))
* remove snakemke module and update snakemake version in requirements ([56f3c05](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/56f3c054fc210b5b662ffe6650537fc98309de96))
* remove unneeded lines ([080a86b](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/080a86bbcf4479ab94c1f56395325af0a9f048c7))
* **Snakefile:** set gvcf_records output and remove unused chr wildcard ([9cabdf6](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/9cabdf6068b3300df7ba3402bbe5a4b2d7d55e2c))
* throw error if no samples found in samplesheet ([30dc3fc](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/30dc3fcea5294177a5915f7dd0125f30606fa341))
* update and rename start_Poirot_marvin.sh to start_HG_marvin.sh ([154e140](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/154e140457464fc9ceac243e6c86611913928c85))
* update cnv_sv module ([1bc9357](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/1bc9357ff59ffc05b4d7934629acdd19dce853c3))
* update coverage file ([9b654a6](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/9b654a69c856b64c6649ac4747f85ae5a1b7c451))
* Update create_excel.py ([ed8bd36](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/ed8bd361afa3d6e056452bb6b045b234e1e72c46))
* Update input files for multiqc report ([84550c3](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/84550c3e14318bb333f4df6654bab689e29dd356))
* update multiqc_config_DNA.yaml ([66cdd5f](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/66cdd5f176f75bae78b3e313f40bd12ece04d730))
* update random_stuff ([6c73a97](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/6c73a9790b351b744d630b7208b65fc239fcc176))
* WGS genepanel in config.yaml ([d1018a9](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/d1018a9811f52fd7069cd6eed232904dffdf9eaa))


### Performance Improvements

* fix sample order so it works for up to 999 samples ([d9344d3](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/d9344d38bd01c96122cc5409e7f53922d4a5f119))


### Documentation

* add a rulegraph for the pipeline ([4f9ac0e](https://www.github.com/clinical-genomics-uppsala/poirot_rd_wgs/commit/4f9ac0e774ca68f41634f9a3f84146c4eed72425))

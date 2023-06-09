# Poirot RD WGS
 Clinical Genomics Uppsala inheritance disease pipeline for WGS made as a snakemake workflow.


The pipeline is built to analys WGS data. Where possible, hydra-genetics modules (https://github.com/hydra-genetics) is being used. The main parts are the same as the GMS nextflow pipeline https://github.com/nf-core/raredisease.


**SNV and indel analysis**

- fastq to BAM with bwa and marking duplicates (https://docs.nvidia.com/clara/parabricks/4.1.0/documentation/tooldocs/man_fq2bam.html#man-fq2bam)
- deepVariant (+ GLNexus for peddy) for calling (https://docs.nvidia.com/clara/parabricks/4.1.0/documentation/tooldocs/man_deepvariant.html#man-deepvariant)


**CNV, and other SV: inversions, deletion and duplications**

- Manta
- CNVpytor
- tiddit
- When possible, we will continue buildning this part, adding other callers maybe CNVkit and delly.
- Combine the results from different callers: SVDB to one vcf-file
  - SVDB will help remove false positives
- Region Of Homozygosity and UniParental Disomy
  - AutoMap (https://github.com/mquinodo/AutoMap) and https://github.com/bjhall/upd
- SMNCopyNumberCaller (https://github.com/Illumina/SMNCopyNumberCaller, https://www.nature.com/articles/s41436-020-0754-0?proof=t)
  - Maybe look into: SMNca (https://onlinelibrary.wiley.com/doi/full/10.1002/humu.24120) and other ways to handle SMN1 och SMN2?

**Mitochondria**

- heteroplasmy (sensitivity) 


**Repeat expansions**

- ExpansionHunter
- annotater with STRanger
- REViewer makes histogram with size distribution per sample


**QC**

- MultiQC report
- Coverage for genes and gene panels
- Kinship and sex-check with peddy


**To implement**

- RNA
- GATK CNV germline caller
- Continued work on SV calling
- Mobile elements
- Several sex-checks
  - samtools idxstats helps with determining sex, can see XXY and females with highly homozygote chrX (make a table with predicted sex based on this)

*Maybe in future*
- Telomerecat is a tool for estimating the average telomere length (TL) for a paired end, whole genome sequencing (WGS) sample
- Cyrius for good call of CYP2D6

# Poirot RD WGS
 Clinical Genomics Uppsala inheritance disease pipeline for WGS made as a snakemake workflow.


The pipeline is built to analys WGS data. Where possible, hydra-genetics modules (https://github.com/hydra-genetics) is being used. The main parts are the same as the GMS nextflow pipeline https://github.com/nf-core/raredisease.


**SNV and indel analysis**

- fastq to BAM with bwa and marking duplicates (https://docs.nvidia.com/clara/parabricks/4.0.0/Documentation/ToolDocs/man_fq2bam.html#man-fq2bam)
- deepVariant (+ GLNexus for peddy) for calling


**CNV, and other SV: inversions, deletion and duplications**

- Manta
- CNVpytor
- tiddit
- When possible, we will continue buildning this part, adding other callers maybe CNVkit and delly.
- Combine the results from different callers: SVDB to one vcf-file
  - SVDB will help remove false positives
- To implement: Region Of Homozygosity and UniParental Disomy
  - AutoMap (https://github.com/mquinodo/AutoMap) and https://github.com/bjhall/upd


**Repeat expansions**

- ExpansionHunter
- annotater with STRanger
- REViewer makes histogram with size distribution per sample


**QC**

- MultiQC report
- To implement: coverage for gene panels
- kinship and sex-check with peddy
 - samtools idxstats helps with determining sex, can see XXY and females with highly homozygote chrX
 

**To implement: SMA**

- SMNCopyNumberCaller (https://github.com/Illumina/SMNCopyNumberCaller, https://www.nature.com/articles/s41436-020-0754-0?proof=t)
- SMNca (https://onlinelibrary.wiley.com/doi/full/10.1002/humu.24120)
- other ways to handle SMN1 och SMN2?


**To implement: Mitochondria**
- heteroplasmy (sensitivity) 


**To implement RNA**


---

#Software or thoughts for future

- **Telomerecat is a tool for estimating the average telomere length (TL) for a paired end, whole genome sequencing (WGS) sample** (Panos kanske är intresserad av svaret)
- Cyrius for good call of CYP2D6

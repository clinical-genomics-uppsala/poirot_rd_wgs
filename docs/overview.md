# Poirot RD WGS

## Mapping reads and BAM file processing

When GPUs are available Poirot can be configured to use [Nvidia's Parabricks](https://www.nvidia.com/en-gb/clara/parabricks/) for read mapping using [fq2bam](https://docs.nvidia.com/clara/parabricks/latest/documentation/tooldocs/man_fq2bam.html#man-fq2bam) tool. This tool performs read mapping with a GPU-accelerated version of BWA-mem, sorting and marking of duplicates.

When only CPUs are available Poirot can be configured perform the read mapping, sorting and duplicate marking on CPU.

- read mapping [BWA-MEM](https://github.com/lh3/bwa)
- read sorting [Samtools sort](https://www.htslib.org/doc/samtools-sort.html)
- marking duplicates with [Picard MarkDuplicates](https://broadinstitute.github.io/picard/command-line-overview.html#MarkDuplicates)

## Variant Calling

### SNV and INDELs

- [Parabricks DeepVariant](https://docs.nvidia.com/clara/parabricks/latest/documentation/tooldocs/man_deepvariant.html#man-deepvariant) when run on GPU or [Google's DeepVariant](https://github.com/google/deepvariant) when run on CPU
- [Glnexus](https://github.com/dnanexus-rnd/GLnexus)
    - Used to create a multisample VCF file analysed with Peddy.
    - Used for the creation of trio VCF files used for UPD analysis

### Mitochondrial short variants

- [GATK Mitochondrial short variant discovery](https://gatk.broadinstitute.org/hc/en-us/articles/4403870837275-Mitochondrial-short-variant-discovery-SNVs-Indels) (SNVs + Indels)

### CNVs and SVs

- CNV callers
    - [CNVpytor](https://github.com/abyzovlab/CNVpytor)

- Structural variant callers
    - [Manta](https://github.com/Illumina/manta)
    - [Tiddit](https://github.com/SciLifeLab/TIDDIT)

- Merging and filtering of SV VCF files
    
    - [SVDB merge](https://github.com/J35P312/SVDB?tab=readme-ov-file#merge) used to merge the Tiddit, Manta and CNVpytor VCF files
    - [SVDB query](https://github.com/J35P312/SVDB?tab=readme-ov-file#query) used to annotate the merge VCF with frequency information from local SV databases
    - Annotation of SVDB merged VCF with Gnomad v4.0 AF using the [Ensembl Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html)
    - Filtering of SV annotated VCF files based on Gnomad AF and the frequency of each SV called in local svdb databases

### Repeat expansions

- Calling and QC of repeat expansions calls
    - [Expansion Hunter](https://github.com/Illumina/ExpansionHunter) estimates the size of short tandem repeats from WGS illumina data.
- QC of repeat expansions calls
    - [REViewer](https://github.com/Illumina/REViewer) is a tools for visualising the read support for the calls made by Expansion Hunter
- Annotation and determination of pathogenic status of repeats with [STRanger](https://github.com/Clinical-Genomics/stranger)

### Regions Of Homozygosity
- [AutoMap](https://github.com/mquinodo/AutoMap)

### SMN Copy Number
 - [SMNCopyNumberCaller](https://github.com/Illumina/SMNCopyNumberCaller)

### UniParental Disomy 
- [upd](https://github.com/bjhall/upd)

## QC

Poirot produces a MultiQC-report for the entire sequencing run to enable easier QC tracking. The report starts with a general statistics table showing the most important QC-values followed by additional QC data and diagrams. The entire MultiQC html-file is interactive and you can filter, highlight, hide or export data using the ToolBox at the right edge of the report.

- The MultiQC-report contains QC data from the following programs:
    - [fastqc](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
    - [mosdepth](https://github.com/brentp/mosdepth)
    - [peddy](https://github.com/brentp/peddy)
    - [Picard CollectAlignmentSummaryMetrics](https://broadinstitute.github.io/picard/command-line-overview.html#CollectAlignmentSummaryMetrics)
    - [Picard CollectDuplicateMetrics](https://gatk.broadinstitute.org/hc/en-us/articles/360042915371-CollectDuplicateMetrics-Picard)
    - [Picard CollectHsMetrics](https://broadinstitute.github.io/picard/command-line-overview.html#CollectHsMetrics)
    - [Picard CollectGcBiasMetrics](https://broadinstitute.github.io/picard/command-line-overview.html#CollectGcBiasMetrics)
    - [Picard CollectInsertSizeMetrics](https://broadinstitute.github.io/picard/command-line-overview.html#CollectInsertSizeMetrics)
    - [Picard CollectWgsMetrics](https://broadinstitute.github.io/picard/command-line-overview.html#CollectWgsMetrics)
    - [samtools stats](https://www.htslib.org/doc/samtools-stats.html)
    - [samtools idxstats](https://www.htslib.org/doc/samtools-idxstats.html)
    - [verifybamid2](https://github.com/Griffan/VerifyBamID)

- Coverage for genes and gene panels. Results written to an excel spreadsheet with a tab for each gene panel.

**To implement**

- GATK CNV germline caller
- Continued work on SV calling and filtering
- Mobile elements
- Several sex-checks
  - samtools idxstats helps with determining sex, can see XXY and females with highly homozygote chrX (make a table with predicted sex based on this)


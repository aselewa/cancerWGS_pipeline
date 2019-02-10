
## introduction

This is a Snakemake pipeline written for the processing of whole-genome sequencing data of matched normal-tumor samples. The pipeline takes as input Illumina FASTQ files and will output:

* Germline SNPs and tumor BAF at these positions (HaplotypeCaller)
* Simple somatic mutations (Mutect2)
* Somatic CNAs (TitanCNA)

The pipeline uses a combination of GATK4 and TitanCNA for calling somatic mutations and copy number alterations.

## Environment

The environment is handled by the Anaconda package manager

#!/bin/python

import glob
import os

pd = config["proj_dir"]
fastq = pd + "fastq/"
aligned = pd + "aligned/"
output =  pd + "output/"
VCFs = output + "VCFs/"

# Make necessary directories
dir_log = pd + "log/"
if not os.path.isdir(dir_log):
    os.mkdir(dir_log)

tmp = pd + "tmp/"
if not os.path.isdir(tmp):
    os.mkdir(tmp)

## Files
REF_GENOME = config["REF_GENOME_DIR"]

## Global wildcards

#patient id
id = glob_wildcards(fastq + "{id}_tumor_R1.fastq.gz").id

#chromosomes
chr = ["%.1d" % i for i in range(1,23)] + ['X','Y']

#types
type = ['tumor','normal']


####### BEGIN
    
rule all:
    input:
        expand(aligned + "{i}_{t}.bam",i=id,t=type),    # ALIGNED BAM FILE
        expand(VCFs + "{i}_germline_het_SNPs.vcf.gz",i=id),            # GERMLINE HETERZYGOUS SNPs
        expand(VCFs + "{i}_somatic_SNPs_filtered.vcf.gz",i=id),        # SOMATIC SNPs
        expand(output + "{i}_tumor_Titan.txt",i=id),                     # PROCESSED TUMOR for TITAN CNV Calling    
        expand(output + "{i}_{t}_readCounts_1k_Titan.wig",i=id,t=type),              # READ COUNTS OF TUMOR/NORMAL for TITAN CNV Calling
        expand(output + "{i}_TitanCNA_Calls.txt",o=output,i=id),
        expand(pd + "{i}_confirm_finish.txt",i=id)

rule bwa:
    input:
        fastq + "{id}_{type}_R1.fastq.gz",
        fastq + "{id}_{type}_R2.fastq.gz"
    output:
        aligned + "{id}_{type}.bam"
    threads: 16
    shell:
        "bwa mem -M -t {threads} {REF_GENOME} {input[0]} {input[1]} | samtools view -bS - > {output}"

rule sort_bams:
    input:
        aligned + "{id}_{type}.bam"
    output:
        output + "{id}_{type}_sorted.bam"
    threads: 10
    shell:
        "samtools sort -@ {threads} -o {output} -O bam {input} -T {output}_temp"

rule remove_duplicates:
    input:
        aligned + "{id}_{type}_sorted.bam"
    output:
        output + "{id}_{type}_sorted_dedupped.bam"
    shell:
        "picard MarkDuplicates I={input} O={output} M={output}_metrics.txt -Xmx8G REMOVE_DUPLICATES=TRUE TMP_DIR={tmp} && samtools index {output}"

# QUALITY BASE RECALIBRATION (OPTIONAL)
#rule base_calibrate:
#    input:
#        output + "{sample}_sorted_dedupped.bam", 
#        output + "{sample}_sorted_dedupped.bam.bai"
#    output:
#        output + "{sample}_recal_data.table"
#    shell:
#        "gatk --java-options '-Xmx8G -Djava.io.tmpdir={tmp}' BaseRecalibrator -I {input[0]} -R {REF_GENOME} -O {output} --known-sites {sites_variation}Mills_and_1000G_gold_standard.indels.b37.vcf --known-sites {sites_variation}10# 00G_phase1.indels.b37.vcf"

#rule recalibrate_bam:
#     input:
#        output + "{sample}_sorted_dedupped.bam", 
#	output + "{sample}_recal_data.table"
#     output:
#        output + "{sample}_sorted_dedupped_calib.bam"
#     shell:
#       "gatk --java-options '-Xmx8G -Djava.io.tmpdir={tmp}' ApplyBQSR -I {input[0]} -R {REF_GENOME} --bqsr-recal-file {input[1]} -O {output}" 


## Call Germline SNPs and prepare inputs for TITAN

rule call_germline_snps:
    input:
        output + "{id}_normal_sorted_dedupped.bam"
    output:
        temp(VCFs + "{id}_{chr}_germline_SNPS_unfiltered.vcf.gz")
    params:
        chrom="{chr}"
    shell:
        "gatk --java-options '-Xmx1G -Djava.io.tmpdir={tmp}' HaplotypeCaller -I {input} -R {REF_GENOME}  --genotyping-mode DISCOVERY -O {output} -L {params.chrom}"

rule concat_filter_vcf_files:
    input:
        expand(VCFs + "{i}_{c}_germline_SNPS_unfiltered.vcf.gz",i=id,c=chr)
    output:
        VCFs + "{id}_germline_het_SNPs.vcf.gz"
    shell:
        "bcftools concat {input} | bcftools view --types snps | bcftools filter --include \"GT='0/1'\" | gzip -c > {output}"

rule get_germline_SNP_pos:
    input:
        VCFs + "{id}_{chr}_germline_SNPS_unfiltered.vcf.gz"
    output:
        temp(VCFs + "{id}_{chr}_germline_SNPs_VARIANT_POS.txt")
    shell:
        "bcftools query -f '%CHROM \t %POS' {input} > {output}"
        
rule tumor_VCF_at_germline_SNPs:
    input:
        VCFs + "{id}_{chr}_germline_SNPs_VARIANT_POS.txt",
        output + "{id}_tumor_sorted_dedupped.bam"
    output:
        temp(VCFs + "{id}_{chr}_tumor_at_germline_SNPs.vcf.gz")
    shell:
        "bcftools mpileup -f {REF_GENOME} -R {input[0]} -o {output} -O z {input[1]}"
    
rule concat_tumor_vcf_files:
    input:
        expand(VCFs + "{i}_{c}_tumor_at_germline_SNPs.vcf.gz",i=id,c=chr)
    output:
        output + "{id}_tumor_Titan.txt"
    shell:
        "bcftools concat {input} | bcftools view --types snps | bcftools query -f '%CHROM,%POS,%REF,%ALT,[%I16\n]' - | awk -F ',' '{{ref=$5+$6; alt=$7+$8; print $1,$2,$3,ref,$4,alt}}' | sed 's/ /\t/g' - > {output}"

rule readCounts_in_window:
    input:
        output + "{id}_{type}_sorted_dedupped.bam"
    output:
        output + "{id}_{type}_readCounts_1k_Titan.wig"
    shell:
        "readCounter -w 1000 {input} > {output}"

rule call_somatic_cnvs:
    input:
        output + "{id}_tumor_Titan.txt",
        output + "{id}_tumor_readCounts_1k_Titan.wig",
        output + "{id}_normal_readCounts_1k_Titan.wig",
    output:
        output + "{id}_TitanCNA_Calls.txt"
    params:
        GC_FILE=config["GC_FILE"],
        MAP_FILE=config["MAP_FILE"],
        nclust=config["CLONAL_CLUSTERS"]
    shell:
        "Rscript scripts/titanCNA_analysis.R {input[0]} {input[1]} {input[2]} {params.GC_FILE} {params.MAP_FILE} {params.nclust} {output}" 

## Call Somatic Mutations

rule call_somatic_mutations:
    input:
        output + "{id}_tumor_sorted_dedupped.bam",
        output + "{id}_normal_sorted_dedupped.bam"
    output:
        temp(VCFs + "{id}_{chr}_somatic_snps_indels.vcf.gz")
    params:
        chrom="{chr}",
        tumor=config["tumor_sample_name"],
        normal=config["normal_sample_name"]
    shell:
        "gatk --java-options '-Xmx1G -Djava.io.tmpdir={tmp}' Mutect2 -I {input[0]} -tumor {params.tumor} -I {input[1]} -normal {params.normal} -R {REF_GENOME} -O {output} -L {params.chrom}"

rule concat_somatic_vcfs:
    input:
        expand(VCFs + "{i}_{c}_somatic_snps_indels.vcf.gz",i=id,c=chr)
    output:
        VCFs + "{id}_somatic_SNPs_Indels_unfiltered.vcf.gz"
    shell:
        "bcftools concat {input} | bgzip -c > {output} && tabix -p vcf {output}"

rule filter_somatic_calls:
    input:
        VCFs + "{id}_somatic_SNPs_Indels_unfiltered.vcf.gz"
    output:
        VCFs + "{id}_somatic_SNPs_Indels_filtered.vcf.gz"
    shell:
        "gatk FilterMutectCalls -V {input} -O {output}"

rule keep_pass_SNPs:
    input:
        VCFs + "{id}_somatic_SNPs_Indels_filtered.vcf.gz"
    output:
        VCFs + "{id}_somatic_SNPs_filtered.vcf.gz"
    params:
        filter="PASS"
    shell:
        "gatk SelectVariants -R {REF_GENOME} --variant {input} -O {output} --select-type-to-exclude INDEL --restrict-alleles-to BIALLELIC --exclude-filtered"

rule clean_up_and_finish:
     input:
        VCFs + "{id}_somatic_SNPs_filtered.vcf.gz",
        VCFs + "{id}_tumor_Titan.txt",
     output:
        pd + "{id}_confirm_finish.txt"
     shell:
        "rm -r {VCFs}*.tbi {VCFs}{wildcards.id}_splitted* && echo 'DONE!' > {output}"
 
#### COMPLETE


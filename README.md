# Chimpanzee Genomics Pilot Project - analysis of whole-genome short-read DNA-seq data

This repository contains scripts for performing QC, preprocessing data, aligning raw reads, variant calling, and preliminary population genetic analysis of WGS data from six West African wild chimpanzees.

## 1- QC

## 2- Mapping

This folder contains codes to prepare mapped BAMs and a small Nextflow pipeline for GATK-based processing.

Prerequisites
- Access to an HPC (for the MPI-PL, lux13 and Gridmaster).
- samtools, bowtie2, gatk, picard, java, bcftools, mosdepth installations (or module availability on HPC).

Scripts and their purposes
- mapping/index_genome.sh  
  Builds the Bowtie2 index for the reference genome (runs `bowtie2-build`). You need to run this once per reference genome.

- mapping/alignment.sh  
  Main read alignment script. Runs `bowtie2` → `samtools view` → `samtools sort` → `samtools index`, creates per-sample job scripts, and submits each job to the Gridmaster using `qsub`. It outputs sorted BAMs and indexes in your "OUTDIR".

- mapping/merge_bams.sh  
  Merges BAM files that belong to the same biological sample (e.g., multiple seq. lanes for the same sample) using `samtools merge`, and indexes the merged BAMs.

- mapping/mark_duplicates.sh  
  Submits jobs that call `gatk MarkDuplicates` on sorted and merged BAMs, producing duplicate-marked BAMs and relevant metrics. You can check the logs/metrics in the `logs/` folder.

- mapping/add_read_groups.sh  
  Adds read-group (RG) information using Picard's `AddOrReplaceReadGroups`. It outputs `.addRG.marked.sorted.bam` files, which are needed for downstream GATK analyses.

- mapping/gatk.nf  
- mapping/test.nf  
  Ignore these Nextflow scripts. I tried to implement the pipeline using Nextflow to automate things, but I need more time to finalise this.

Run order:
1. mapping/index_genome.sh — build Bowtie2 index (if missing).
2. mapping/alignment.sh — align raw reads → sorted BAMs.
3. mapping/merge_bams.sh — merge lane BAMs per sample (if needed).
4. mapping/mark_duplicates.sh — mark duplicates on sorted/merged BAMs.
5. mapping/add_read_groups.sh — add read groups for GATK.

Note: Make sure to check and edit inDir/outDir/ref etc. variables before running each script, as most directories are hard-coded.
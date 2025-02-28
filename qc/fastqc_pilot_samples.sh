#!/bin/bash

# This script will run fastqc on raw fastq files.
# Gokberk Alagoz, 28.02.25
##########################

# Paths
inDir="/data/workspaces/lag/workspaces/lg-evolution/working_data/tai_chimp_genomics/tai_pilot_data/X204SC24064024-Z01-F001/01.RawData"
outDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/fastqc/pilot_samples"

# Run fastqc
module load fastqc/0.11.9

fastqc -o ${outDir} ${inDir}/*/*fq.gz

#!/bin/bash

# This script will merge bam files that belong to
# same biological samples. Sometimes, the sequencer
# generates multiple fastq files for the same sample.
# This may be because of lane or run splitting, which
# sequencing companies do allow for more samples to be
# sequenced in a single run.
#
# Gokberk Alagoz - 11.03.25
###########################

# Paths and modules
inDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/mapping/"

module purge
module load samtools/1.9

# Merge bams
# TAI3468
samtools merge ${inDir}TAI3468_MKDN240009949-1A_E250024699_E250025014_L1.sorted.bam ${inDir}TAI3468_MKDN240009949-1A_E250024699_L1.sorted.bam ${inDir}TAI3468_MKDN240009949-1A_E250025014_L1.sorted.bam
samtools index ${inDir}TAI3468_MKDN240009949-1A_E250024699_E250025014_L1.sorted.bam

# TAI661
samtools merge ${inDir}TAI661_MKDN240009947-1A_E200024585_E250025003_L1.sorted.bam ${inDir}TAI661_MKDN240009947-1A_E200024585_L1.sorted.bam ${inDir}TAI661_MKDN240009947-1A_E250025003_L1.sorted.bam
samtools index ${inDir}TAI661_MKDN240009947-1A_E200024585_E250025003_L1.sorted.bam


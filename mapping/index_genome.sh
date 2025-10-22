#!/bin/bash
# This script will index the chimpanzee genome.
# GA - 06.03.25

# Paths

inDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/resources/genomes/chimpanzee/NHGRI_mPanTro3-v2.0_pri_ncbi/ncbi_dataset/data/GCF_028858775.2/"

# Index

bowtie2-build ${inDir}GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna panTro3-v2.0

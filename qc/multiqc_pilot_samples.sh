#!/bin/bash

# This script will run multiqc on fastqc results.
# Multiqc does not run any analysis, it makes a nice
# report of fastqc results from multiple samples.
#
# When you run multiqc for the first time, follow the
# manual here to install multiqc and set up a conda 
# environment.
#
# Gokberk Alagoz, 28.02.25
##########################

# Paths
inDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/fastqc/pilot_samples"

# Run multiqc

module load miniconda/3.2021.10
conda activate myenv

cd ${inDir}
multiqc ${inDir} # automatically creates a file named multiqc_report.html in the current working directory

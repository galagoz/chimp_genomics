#!/bin/bash
# This script will calculate the sequencing depth,
# and will allow variant filtering based on positional
# depth. Following Kuderna et al. (2023), I use mosdepth
# package for this: https://github.com/brentp/mosdepth
# Gokberk Alagoz - 10.04.2025

# Paths

inDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/haplotype_calling"
outDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/mapping/sequencing_depths"

# Make directories
mkdir -p "$outDir"
mkdir -p "${outDir}/logs"

# Find all BAM files
find "$inDir" -type f -name "*.addRG.marked.REsorted.bam" | while read -r BAM_FILE; do
    # Extract sample paramaeters
    SAMPLE_ID=$(basename "$BAM_FILE" | sed 's/.addRG.marked.REsorted.bam//')

    # Create a job script for each sample
    JOB_SCRIPT="${outDir}/logs/mosdepth_job_${SAMPLE_ID}.sh"

    cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#$ -N seqDepthCalc_${SAMPLE_ID}
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash
#$ -o ${outDir}/logs/mosdepth_${SAMPLE_ID}.out     # STDOUT log
#$ -e ${outDir}/logs/mosdepth_${SAMPLE_ID}.err     # STDERR log

echo "Processing sample: $SAMPLE_ID"

module load miniconda/3.2021.10
conda activate mosdepthEnv

mosdepth --fast-mode --by 500 ${SAMPLE_ID} ${BAM_FILE}

# --fast-mode avoids the extra calculations of mate pair overlap and cigar operations,
#and also allows htslib to extract less data from CRAM, providing a substantial speed improvement.

echo "Finished processing $SAMPLE_ID"
EOF

    # Submit the job to the grid
    qsub "$JOB_SCRIPT"

    echo "Submitted job for $SAMPLE_ID"
done

echo "All jobs submitted!"


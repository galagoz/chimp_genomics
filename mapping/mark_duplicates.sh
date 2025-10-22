#!/bin/bash
# This script will mark duplicate reads in bam files.
# This is needed in case the duplicate artifacts are
# incorrectly detected as multiple clusters by the
# sequencer. (see GATK documentation for further info)
# Gokberk Alagoz - 17.03.2025

# Paths
READS_DIR="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/mapping"
OUTDIR="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/mapping"

# Make directories
mkdir -p "$OUTDIR"
mkdir -p "${OUTDIR}/logs"

# Find all BAM files
find "$OUTDIR" -type f -name "*.sorted.bam" | while read -r BAM_FILE; do
    # Extract sample ID
    SAMPLE_ID=$(basename "$BAM_FILE" | sed 's/.sorted.bam//')

    # Define output files
    MARKED_BAM="${OUTDIR}/${SAMPLE_ID}.marked.sorted.bam"
    METRICS="${OUTDIR}/${SAMPLE_ID}.metrics.txt"

    # Create a job script for each sample
    JOB_SCRIPT="${OUTDIR}/logs/job_${SAMPLE_ID}.sh"

    cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#$ -N mark_duplicates_${SAMPLE_ID}        # Job name
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash
#$ -o ${OUTDIR}/logs/${SAMPLE_ID}.out     # STDOUT log
#$ -e ${OUTDIR}/logs/${SAMPLE_ID}.err     # STDERR log

module load java/jdk1.8.0_201 gatk/4.3.0.0

echo "Processing sample: $SAMPLE_ID"

gatk MarkDuplicates -I "$BAM_FILE" -O "$MARKED_BAM" -M "$METRICS"

echo "Finished processing $SAMPLE_ID"
EOF

    # Submit the job to the grid
    qsub "$JOB_SCRIPT"

    echo "Submitted job for $SAMPLE_ID"
done

echo "All jobs submitted!"

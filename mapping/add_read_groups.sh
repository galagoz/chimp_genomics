#!/bin/bash
# This script will add reads groups to bam files.
# Each read should have a read group so that variant calling
# can be performed properly.
# (see GATK documentation for further info)
# Gokberk Alagoz - 17.03.2025

# Paths
inDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/mapping"
outDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/mapping/read_groups"

# Make directories
mkdir -p "$outDir"
mkdir -p "${outDir}/logs"

# Find all sorted and duplicate marked BAM files
find "$inDir" -type f -name "*.marked.sorted.bam" | while read -r BAM_FILE; do
    # Extract sample paramaeters
    SAMPLE_ID=$(basename "$BAM_FILE" | sed 's/.marked.sorted.bam//')
    SAMPLE_RGID=$(basename "$BAM_FILE" | sed -r 's/.*1A_(.*?).marked.*/\1/')
    LIBRARY_ID=$(basename "$BAM_FILE" | cut -d '_' -f 2)
    ANIMAL_ID=$(basename "$BAM_FILE" | cut -d '_' -f 1)

    # Define output files
    RG_BAM="${outDir}/${SAMPLE_ID}.addRG.marked.sorted.bam"

    # Create a job script for each sample
    JOB_SCRIPT="${outDir}/logs/job_${SAMPLE_ID}.sh"

    cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#$ -N add_read_groups_${SAMPLE_ID}
#$ -cwd
#$ -q single15.q
#$ -S /bin/bash
#$ -o ${outDir}/logs/${SAMPLE_ID}.out     # STDOUT log
#$ -e ${outDir}/logs/${SAMPLE_ID}.err     # STDERR log

module load java/openjdk-17.0.2 picard/3.3.0

echo "Processing sample: $SAMPLE_ID"

picard AddOrReplaceReadGroups -I ${BAM_FILE} \
			      -O ${RG_BAM} \
			      -RGID ${SAMPLE_RGID} \
			      -RGLB ${LIBRARY_ID} \
			      -RGPL DNBSEQ \
			      -RGPU ${SAMPLE_RGID}_${LIBRARY_ID} \
			      -RGSM ${ANIMAL_ID}

echo "Finished processing $SAMPLE_ID"
EOF

    # Submit the job to the grid
    qsub "$JOB_SCRIPT"

    echo "Submitted job for $SAMPLE_ID"
done

echo "All jobs submitted!"

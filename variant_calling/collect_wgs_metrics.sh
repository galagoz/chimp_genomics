#!/bin/bash
# This script will collect WGS metrics about coverage
# and performance of WGS experiments from BAM files.
# Gokberk Alagoz - 29.04.2025

# Paths

inDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/mapping"
refFile="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/resources/genomes/chimpanzee/NHGRI_mPanTro3-v2.0_pri_ncbi/ncbi_dataset/data/GCF_028858775.2/GCF_028858775.2_NHGRI_mPanTro3-v2.0_pri_genomic.fna"
outDir="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/mapping/WGS_metrics"

# Make directories
mkdir -p "$outDir"
mkdir -p "${outDir}/logs"

# Find all BAM files
find "$inDir" -type f -name "*_L1.sorted.bam" | while read -r BAM_FILE; do
    # Extract sample paramaeters
    SAMPLE_ID=$(basename "$BAM_FILE" | sed 's/._L1.sorted.bam//')
    LIBRARY_ID=$(basename "$BAM_FILE" | cut -d '_' -f 2)
    ANIMAL_ID=$(basename "$BAM_FILE" | cut -d '_' -f 1)

    # Define output files
    WGS_METRICS="${outDir}/${SAMPLE_ID}_wgs_metics.txt"

    # Create a job script for each sample
    JOB_SCRIPT="${outDir}/logs/job_${SAMPLE_ID}.sh"

    cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#$ -N collectWGSmetrics_${SAMPLE_ID}
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash
#$ -o ${outDir}/logs/${SAMPLE_ID}.out     # STDOUT log
#$ -e ${outDir}/logs/${SAMPLE_ID}.err     # STDERR log

module load java/jre1.8.0_201 gatk/4.3.0.0

echo "Processing sample: $SAMPLE_ID"

gatk CollectWgsMetrics -I ${BAM_FILE} \
		       -R ${refFile} \
                       -O ${WGS_METRICS}

echo "Finished processing $SAMPLE_ID"
EOF

    # Submit the job to the grid
    qsub "$JOB_SCRIPT"

    echo "Submitted job for $SAMPLE_ID"
done

echo "All jobs submitted!"

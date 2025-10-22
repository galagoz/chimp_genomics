#!/bin/bash

# This script will align raw DNA-seq reads
# to the chimpanzee reference genome.
# Gokberk Alagoz - 06.03.25

# Paths
READS_DIR="/data/workspaces/lag/workspaces/lg-evolution/working_data/tai_chimp_genomics/tai_pilot_data/X204SC24064024-Z01-F001/01.RawData"
GENOME_INDEX="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/resources/genomes/chimpanzee/NHGRI_mPanTro3-v2.0_pri_ncbi/ncbi_dataset/data/GCF_028858775.2/panTro3-v2.0"
OUTDIR="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/mapping"

# Make directories
mkdir -p "$OUTDIR"
mkdir -p "logs"

# Find all forward read files (_1.fq.gz)
find "$READS_DIR" -type f -name "*_1.fq.gz" | while read -r R1; do
    # Extract sample ID
    SAMPLE_ID=$(basename "$R1" | sed 's/_1.fq.gz//')

    # Define reverse read file
    R2="${R1/_1.fq.gz/_2.fq.gz}"

    # Check if reverse read exists
    if [[ ! -f "$R2" ]]; then
        echo "WARNING: Paired file not found for $R1, skipping..."
        continue
    fi

    # Output BAM file
    BAM_FILE="${OUTDIR}/${SAMPLE_ID}.sorted.bam"

    # Create a job script for this sample
    JOB_SCRIPT="logs/job_${SAMPLE_ID}.sh"

    cat <<EOF > "$JOB_SCRIPT"
#!/bin/bash
#$ -N align_${SAMPLE_ID}        # Job name
#$ -cwd
#$ -q multi15.q
#$ -S /bin/bash
#$ -o logs/${SAMPLE_ID}.out     # STDOUT log
#$ -e logs/${SAMPLE_ID}.err     # STDERR log

module load libtbb/2019 bowtie2/2.4.2 samtools/1.9

echo "Processing sample: $SAMPLE_ID"

bowtie2 -x "$GENOME_INDEX" -1 "$R1" -2 "$R2" | \
samtools view -bS - | \
samtools sort -o "$BAM_FILE"

samtools index "$BAM_FILE"

echo "Finished processing $SAMPLE_ID -> $BAM_FILE"
EOF

    # Submit the job
    qsub "$JOB_SCRIPT"

    echo "Submitted job for $SAMPLE_ID"
done

echo "All jobs submitted!"


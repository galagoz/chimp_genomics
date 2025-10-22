#!/bin/bash
# This script will calculate heterozygosity for each sample
# and create a summary file with info from all samples.
# Gokberk Alagoz - 23.04.2025

# Set paths
IN_DIR="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/filtered_variants"
CALL_MASK_DIR="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/callability_masks"
OUT_DIR="/data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/heterozygosity_rates"
LOG_DIR="${OUT_DIR}/logs"
SUMMARY_FILE="${OUT_DIR}/heterozygosity_summary.tsv"

# Create output dirs
mkdir -p "$OUT_DIR" "$LOG_DIR"

# Samples to process
samples=(
  "TAI3468"
  "TAI48406"
  "TAI48807"
  "TAI48879"
  "TAI661"
)

# Create summary file header
echo -e "Sample\tHet_Count\tCallable_Sites\tGenome_Wide_Heterozygosity" > "$SUMMARY_FILE"

# Submit one job per sample
for SAMPLE in "${samples[@]}"; do
  INPUT_VCF="${IN_DIR}/${SAMPLE}.filtered.snps.vcf.gz"
  CALL_MASK_VCF="${CALL_MASK_DIR}/${SAMPLE}.callability_mask.vcf.gz"
  TEMP_OUT="${OUT_DIR}/${SAMPLE}.heterozygosity.tmp"

  qsub -N hetero_${SAMPLE} \
       -o "${LOG_DIR}/${SAMPLE}.o" -e "${LOG_DIR}/${SAMPLE}.e" <<EOF
#!/bin/bash
module load bcftools/1.18

SAMPLE=$SAMPLE
INPUT_VCF=$INPUT_VCF
CALL_MASK_VCF=$CALL_MASK_VCF
OUT_DIR=$OUT_DIR
TEMP_OUT=$TEMP_OUT

echo "Processing sample \$SAMPLE"

bcftools filter -e "GT!='het' | (FORMAT/AD[0:1] / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) < 0.25) | (FORMAT/AD[0:1] / (FORMAT/AD[0:0] + FORMAT/AD[0:1]) > 0.75)" \\
    -Oz -o "\$OUT_DIR/\$SAMPLE.filtered_hetAB.vcf.gz" "\$INPUT_VCF" || exit 1

bcftools index -c "\$OUT_DIR/\$SAMPLE.filtered_hetAB.vcf.gz" || exit 1

HET_COUNT=\$(bcftools view -H -i 'GT="het"' "\$OUT_DIR/\$SAMPLE.filtered_hetAB.vcf.gz" | wc -l)
CALLABLE_SITES=\$(bcftools view -H "\$CALL_MASK_VCF" | wc -l)

if [[ "\$CALLABLE_SITES" -gt 0 ]]; then
  HET_RATE=\$(echo "scale=10; \$HET_COUNT / \$CALLABLE_SITES" | bc)
else
  HET_RATE="NA"
fi

echo -e "\$SAMPLE\t\$HET_COUNT\t\$CALLABLE_SITES\t\$HET_RATE" > "\$TEMP_OUT"

echo "Sample \$SAMPLE done."
EOF

done

echo "Submitted all jobs. Waiting for completion..."

# Wait for all jobs to finish
while true; do
  RUNNING=$(qstat -u $USER | grep -c hetero_)
  if [[ "$RUNNING" -eq 0 ]]; then
    break
  fi
  echo "$RUNNING jobs still running..."
  sleep 30
done

echo "All jobs finished. Merging results..."

# Merge all temp results
cat "$SUMMARY_FILE" > "${SUMMARY_FILE}.tmp"
for SAMPLE in "${samples[@]}"; do
  cat "${OUT_DIR}/${SAMPLE}.heterozygosity.tmp" >> "${SUMMARY_FILE}.tmp"
done
mv "${SUMMARY_FILE}.tmp" "$SUMMARY_FILE"

echo "Heterozygosity summary saved to $SUMMARY_FILE"

#!/bin/bash
# This script will calculate heterozygosity.
# Gokberk Alagoz - 23.04.2025 (recevied a lot of help from GPT-4).

# Define samples list, each sample line:
# sample_name input_vcf callability_mask_vcf

samples=(
  "TAI3468 /data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/filtered_variants/TAI3468.filtered.snps.vcf.gz /data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/callability_masks/TAI3468.callability_mask.vcf.gz"
  "TAI48406 /data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/filtered_variants/TAI48406.filtered.snps.vcf.gz /data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/callability_masks/TAI48406.callability_mask.vcf.gz"
  "TAI48807 /data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/filtered_variants/TAI48807.filtered.snps.vcf.gz /data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/callability_masks/TAI48807.callability_mask.vcf.gz"
  "TAI48879 /data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/filtered_variants/TAI48879.filtered.snps.vcf.gz /data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/callability_masks/TAI48879.callability_mask.vcf.gz"
  "TAI661 /data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/filtered_variants/TAI661.filtered.snps.vcf.gz /data/workspaces/lag/workspaces/lg-evolution/analysis/tai_chimp_genomics/results/variant_calling/callability_masks/TAI661.callability_mask.vcf.gz"
)

# Output summary file
OUTFILE="heterozygosity_summary.tsv"

echo -e "Sample\tHet_Count\tCallable_Sites\tGenome_Wide_Heterozygosity" > "$OUTFILE"

module load bcftools/1.18  # adjust as needed for your environment

for line in "${samples[@]}"; do
  read -r SAMPLE INPUT_VCF CALL_MASK_VCF <<< "$line"
  echo "Processing sample: $SAMPLE"

  # Filter heterozygous sites by allele balance (AB = AD_alt / (AD_ref + AD_alt))
  # Exclude heterozygotes with AB < 0.25 or > 0.75
  bcftools filter -e "GT!='het' | (FMT/AD[0:1] / (FMT/AD[0:0] + FMT/AD[0:1]) < 0.25) | (FMT/AD[0:1] / (FMT/AD[0:0] + FMT/AD[0:1]) > 0.75)" \
    -Oz -o "${SAMPLE}.filtered_hetAB.vcf.gz" "$INPUT_VCF"
  bcftools index -c "${SAMPLE}.filtered_hetAB.vcf.gz"

  # Count heterozygous sites passing the allele balance filter
  HET_COUNT=$(bcftools view -H -i 'GT[0]="het"' "${SAMPLE}.filtered_hetAB.vcf.gz" | wc -l)
  echo "Heterozygous count: $HET_COUNT"
  # Count callable sites from callability mask VCF (number of sites)
  # Consider all sites in callability mask VCF as callable
  CALLABLE_SITES=$(bcftools view -H "$CALL_MASK_VCF" | wc -l)
  
  # Calculate genome-wide heterozygosity (avoid division by zero)
  if [ "$CALLABLE_SITES" -gt 0 ]; then
    HET_RATE=$(echo "scale=10; $HET_COUNT / $CALLABLE_SITES" | bc)
  else
    HET_RATE="NA"
  fi

  # Append to output file
  echo -e "${SAMPLE}\t${HET_COUNT}\t${CALLABLE_SITES}\t${HET_RATE}" >> "$OUTFILE"

  # Optional: remove intermediate filtered file or keep for inspection
  # rm -f "${SAMPLE}.filtered_hetAB.vcf.gz" "${SAMPLE}.filtered_hetAB.vcf.gz.csi"
done

echo "All done! Summary saved to $OUTFILE"

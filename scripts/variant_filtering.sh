#!/bin/bash

# Description:
# This script filters somatic variants based on allele frequency thresholds from multiple population datasets
# and extracts variants common between two VCFs (e.g., matched tumor-normal pairs).

# Input files
VCF1="Final.vcf"
VCF2="Test.vcf"

# Output files
SORTED_VCF1="Sorted_Final.vcf"
SORTED_VCF2="Sorted_Test.vcf"
MATCHED_VCF="Matched_Final.vcf"
FILTERED_VCF="Filtered_Matched_Final.vcf"

# Step 1: Sort input VCFs by chromosome (column 1)
sort -k1,1 $VCF1 > $SORTED_VCF1
sort -k1,1 $VCF2 > $SORTED_VCF2

# Step 2: Join matching rows between the two VCFs (assuming matching on chromosome)
join -t $'\t' -1 1 -2 1 $SORTED_VCF1 $SORTED_VCF2 > $MATCHED_VCF

# Step 3: Filter variants with MAF > 0.01 in dbSNP, 1000 Genomes, gnomAD, or other cohorts
awk -F '\t' '
BEGIN {OFS = "\t"}
{
    dbSNP_maf = 0; genome_maf = 0; gnomAD_maf = 0; other_maf = 0;

    for (i = 1; i <= NF; i++) {
        if ($i ~ /dbSNP_MAF=([0-9.]+)/)    dbSNP_maf  = substr($i, RSTART + 10, RLENGTH)
        if ($i ~ /1000genome_MAF=([0-9.]+)/) genome_maf = substr($i, RSTART + 14, RLENGTH)
        if ($i ~ /gnomAD_MAF=([0-9.]+)/)    gnomAD_maf = substr($i, RSTART + 11, RLENGTH)
        if ($i ~ /OtherCohort_MAF=([0-9.]+)/) other_maf = substr($i, RSTART + 15, RLENGTH)
    }

    if (dbSNP_maf > 0.01 || genome_maf > 0.01 || gnomAD_maf > 0.01 || other_maf > 0.01)
        print $0
}' $MATCHED_VCF > $FILTERED_VCF

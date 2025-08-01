##
#!/bin/bash

# ================================
# Title: preprocessing.sh
# Purpose: Raw WGS/WES FASTQ preprocessing for variant calling
# Steps:
#   1. Adapter trimming (Cutadapt)
#   2. BWA alignment (SAM/BAM)
#   3. Mark duplicates & BQSR
#   4. Variant calling (Mutect2) & annotation (Funcotator)
# Author: [Junho Kang]
# ================================

# ========== USER INPUT ==========
TRIM_DIR="./Trimming"
SAM_DIR="./sam"
FILTER_DIR="./filter"
RESOURCE="/path/to/funcotator_dataSources.v1.7.20200521s/"
REF_DIR="/path/to/gatk_bundle/"
GERMLINE="${REF_DIR}af-only-gnomad.hg38.vcf.gz"
REFERENCE="${REF_DIR}Homo_sapiens_assembly38.fasta"

# ========== ENV SETUP ==========
mkdir -p "$TRIM_DIR" "$SAM_DIR" "$FILTER_DIR"

# ========== 1. TRIMMING ==========
echo "[STEP 1] Trimming adapter sequences with Cutadapt..."
for R1 in *_R1.fastq.gz; do
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    cutadapt \
        -a AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
        -A AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT \
        -o "${TRIM_DIR}/${SAMPLE}_R1.fastq.gz" \
        -p "${TRIM_DIR}/${SAMPLE}_R2.fastq.gz" \
        -j 8 \
        "${SAMPLE}_R1.fastq.gz" "${SAMPLE}_R2.fastq.gz"
done

# ========== 2. ALIGNMENT ==========
echo "[STEP 2] Aligning reads using BWA MEM..."
for R1 in ${TRIM_DIR}/*_R1.fastq.gz; do
    SAMPLE=$(basename "$R1" _R1.fastq.gz)
    R2="${TRIM_DIR}/${SAMPLE}_R2.fastq.gz"

    # Run BWA
    bwa mem -M -t 8 -R "@RG\tID:${SAMPLE}\tSM:${SAMPLE}\tLB:${SAMPLE}\tPL:ILLUMINA" \
        "$REFERENCE" "$R1" "$R2" > "${SAM_DIR}/${SAMPLE}.sam"

    # Convert SAM to BAM, sort, and index
    samtools view -bS "${SAM_DIR}/${SAMPLE}.sam" | samtools sort -o "${SAM_DIR}/${SAMPLE}.sorted.bam"
    samtools index "${SAM_DIR}/${SAMPLE}.sorted.bam"
    rm "${SAM_DIR}/${SAMPLE}.sam"
done

# ========== 3. MARK DUPLICATES & BQSR ==========
echo "[STEP 3] Marking duplicates and recalibrating base quality scores..."
for BAM in ${SAM_DIR}/*.sorted.bam; do
    SAMPLE=$(basename "$BAM" .sorted.bam)

    # Mark duplicates
    gatk MarkDuplicatesSpark \
        -I "$BAM" \
        -O "${SAM_DIR}/${SAMPLE}.duplicates.bam" \
        --duplicate-tagging-policy OpticalOnly

    # BQSR steps
    gatk BaseRecalibrator \
        -I "${SAM_DIR}/${SAMPLE}.duplicates.bam" \
        -R "$REFERENCE" \
        --known-sites "${REF_DIR}dbsnp_146.hg38.vcf.gz" \
        --known-sites "${REF_DIR}hapmap_3.3_grch38_pop_stratified_af.vcf.gz" \
        --known-sites "${REF_DIR}Mills_and_1000G_gold_standard.indels.hg38.vcf.gz" \
        -O "${SAM_DIR}/${SAMPLE}.recal_data.table"

    gatk ApplyBQSR \
        -R "$REFERENCE" \
        -I "${SAM_DIR}/${SAMPLE}.duplicates.bam" \
        --bqsr-recal-file "${SAM_DIR}/${SAMPLE}.recal_data.table" \
        -O "${SAM_DIR}/${SAMPLE}.final.bam"
done

# ========== 4. VARIANT CALLING ==========
echo "[STEP 4] Running Mutect2 for somatic variant calling and annotation..."
for BAM in ${SAM_DIR}/*.final.bam; do
    SAMPLE=$(basename "$BAM" .final.bam)

    # Mutect2 + filtering
    gatk Mutect2 \
        -R "$REFERENCE" \
        -I "$BAM" \
        -tumor "$SAMPLE" \
        --germline-resource "$GERMLINE" \
        -O "${SAM_DIR}/${SAMPLE}.vcf.gz"

    gatk FilterMutectCalls \
        -R "$REFERENCE" \
        -V "${SAM_DIR}/${SAMPLE}.vcf.gz" \
        -O "${SAM_DIR}/${SAMPLE}.filtered.vcf.gz"

    # Annotation using Funcotator
    gatk Funcotator \
        -R "$REFERENCE" \
        -V "${SAM_DIR}/${SAMPLE}.filtered.vcf.gz" \
        -O "${SAM_DIR}/${SAMPLE}.maf.gz" \
        --output-file-format MAF \
        --data-sources-path "$RESOURCE" \
        --ref-version hg38 \
        --annotation-default tumor_barcode:$SAMPLE

    # Optional: extract key columns (customize as needed)
    zcat "${SAM_DIR}/${SAMPLE}.maf.gz" | \
        awk -F '\t' 'BEGIN {OFS="\t"} {print $1,$5,$6,$7,$8,$9,$10,$11,$12,$13,$14,$16,$42,$80,$215}' \
        > "${FILTER_DIR}/${SAMPLE}.filtered.tsv"
done

echo "[DONE] All preprocessing steps completed."

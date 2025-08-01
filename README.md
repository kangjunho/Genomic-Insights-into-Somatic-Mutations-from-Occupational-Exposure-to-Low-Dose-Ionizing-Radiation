# Genomic-Insights-into-Somatic-Mutations-from-Occupational-Exposure-to-Low-Dose-Ionizing-Radiation

# Somatic Mutation Analysis Pipeline

This repository contains code for preprocessing, filtering, and visualization of somatic variants derived from WGS data.

## Folder Structure
- `preprocessing.sh` : raw FASTQ trimming and alignment
- `variant_filtering.sh` : post-VCF filtering using MAF
- `visualization.R` : all figures shown in the manuscript (Figure 2 and 3)

## Requirements
- GATK v4
- samtools
- Cutadapt
- R (ggplot2, dplyr, ggpubr, ggbreak, patchwork, tidyr)

## How to Run
1. `bash preprocessing.sh`
2. `bash variant_filtering.sh`
3. `Rscript visualization.R`

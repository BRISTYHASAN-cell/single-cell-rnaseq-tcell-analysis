#!/bin/bash
#SBATCH --job-name=run_all_tcr
#SBATCH --output=out_%j.txt
#SBATCH --error=err_%j.txt
#SBATCH --time=48:00:00
#SBATCH --mem=1550G
#SBATCH --cpus-per-task=4
#SBATCH --partition=amd-hm
#SBATCH --account=pi-rroychowdhury

# Load modules and activate your environment
module load gcc/13.2.0
module load R/4.4.1
module load python/anaconda-2023.09
source activate /scratch/midway3/abmshamimuh/envs_/gen_/

# Set base path
BASE_DIR=/scratch/midway3/abmshamimuh/Stringent
PIPELINE=${BASE_DIR}/vdjfasta-tcr/bin/single-cell-pipeline.sh

# Output directory
OUTPUT_DIR=${BASE_DIR}/TCR_CSVs
mkdir -p ${OUTPUT_DIR}

cd ${BASE_DIR}

# Loop through all R1 FASTQ files
for R1 in *_R1_001.fastq; do
    SAMPLE=$(echo $R1 | cut -d'_' -f1)
    R2="${SAMPLE}_R2_001.fastq"

    if [ -f "$R2" ]; then
        echo "Processing sample: $SAMPLE"
        LOG="${OUTPUT_DIR}/${SAMPLE}_log.txt"
        CSV="${OUTPUT_DIR}/${SAMPLE}_TCR.csv"

        # Run the pipeline
        ${PIPELINE} "$R1" "$R2" > "$LOG"

        # Copy final output to TCR_CSVs
        FINAL=$(find . -type f -name "*.csv" | grep -m 1 "$SAMPLE")
        if [ -f "$FINAL" ]; then
            cp "$FINAL" "$CSV"
        fi
    else
        echo "Missing R2 file for $SAMPLE, skipping..."
    fi
done

# Merge all sample CSVs into one
echo "Merging all sample CSVs..."
head -n 1 $(ls ${OUTPUT_DIR}/*_TCR.csv | head -n 1) > ${BASE_DIR}/all_samples_TCR.csv
for file in ${OUTPUT_DIR}/*_TCR.csv; do
    tail -n +2 "$file" >> ${BASE_DIR}/all_samples_TCR.csv
done

conda deactivate

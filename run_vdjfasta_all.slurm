#!/bin/bash
#SBATCH --job-name=vdjfasta_tcr_all
#SBATCH --output=out_%j.txt
#SBATCH --error=err_%j.txt
#SBATCH --time=24:00:00
#SBATCH --mem=1550G
#SBATCH --cpus-per-task=16
#SBATCH --partition=amd-hm
#SBATCH --account=pi-rroychowdhury

# Load necessary modules
module load gcc/13.2.0
module load R/4.4.1
module load python/anaconda-2023.09

# Activate conda environment
source activate /scratch/midway3/abmshamimuh/envs_/gen_/

# Navigate to the vdjfasta-tcr tool directory
cd /scratch/midway3/abmshamimuh/Stringent/vdjfasta-tcr

# Run the TCR annotation tool
./vdjfasta-tcr.pl \
    -f ../joined.fa \
    -o all_TCR_results.csv \
    -d ../db \
    -m full \
    -id all_samples

# Deactivate environment
conda deactivate

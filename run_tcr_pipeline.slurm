#!/bin/bash
#SBATCH --job-name=annotate_emmaTCR
#SBATCH --output=annotate_emmaTCR_%j.out
#SBATCH --error=annotate_emmaTCR_%j.err
#SBATCH --partition=amd-hm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=01:00:00

echo "Job started at: $(date)"
echo "Activating conda environment..."

# Load Anaconda and activate your environment
module load python/anaconda-2023.09
source activate /scratch/midway3/abmshamimuh/envs_/gen_/

echo "Running vdjfasta-tcr on emma.joined.fa..."
vdjfasta-tcr emma.joined.fa > emmaTCR.csv

echo "Annotation complete."
echo "Output: emmaTCR.csv"
echo "Job finished at: $(date)"

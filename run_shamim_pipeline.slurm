#!/bin/bash
#SBATCH --job-name=shamim_pipeline
#SBATCH --output=out_%j.txt
#SBATCH --error=err_%j.txt
#SBATCH --time=15:00:00
#SBATCH --mem=250G
#SBATCH --cpus-per-task=4
#SBATCH --partition=bigmem
#SBATCH --account=pi-rroychowdhury
 
# Load R module
module load python/anaconda-2023.09
module load R/4.4.1
module load gcc/13.2.0
source activate /scratch/midway3/abmshamimuh/envs_/gen_/
 
# Run R script
Rscript ~/Shamim_Roshni_pipeline.R
 
conda deactivate

#!/bin/bash
#SBATCH --job-name=demux_fasta
#SBATCH --output=out_%j.txt
#SBATCH --error=err_%j.txt
#SBATCH --time=16:00:00
#SBATCH --mem=1550G
#SBATCH --cpus-per-task=16
#SBATCH --partition=amd-hm
#SBATCH --account=pi-rroychowdhury

# Load necessary modules
module load gcc/13.2.0
module load R/4.4.1
module load python/anaconda-2023.09
source activate /scratch/midway3/abmshamimuh/envs_/gen_/

# Move to working directory
cd /scratch/midway3/abmshamimuh/Stringent

# Run demux script
/scratch/midway3/abmshamimuh/vdjfasta-tcr/bin/fasta-demux.pl \
  --dnafile=joined.fa

conda deactivate

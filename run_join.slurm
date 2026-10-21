#!/bin/bash
#SBATCH --job-name=join_fasta
#SBATCH --output=out_%j.txt
#SBATCH --error=err_%j.txt
#SBATCH --time=12:00:00
#SBATCH --mem=100G
#SBATCH --cpus-per-task=4
#SBATCH --partition=bigmem
#SBATCH --account=pi-rroychowdhury

module load gcc/13.2.0
module load R/4.4.1
module load python/anaconda-2023.09
source activate /scratch/midway3/abmshamimuh/envs_/gen_/

# Run join script
/scratch/midway3/abmshamimuh/vdjfasta-tcr/bin/fasta-hiseq-join-slow.pl \
  --forward=cleaned_R1.fasta \
  --reverse=cleaned_R2.fasta \
  --outputfile=joined.fa \
  --ignoredfile=ignored.fa

conda deactivate

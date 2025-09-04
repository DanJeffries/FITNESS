#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=00:30:00 
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=128G
#SBATCH --export=NONE
#SBATCH --array=5
#SBATCH --job-name=SNP_selector_PopSpec
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load VCFtools/0.1.16-GCC-10.3.0

WD=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/analyses/cohort_1/freq/
POPS=FG,LG,SL,SR,TL,WB,WK,WT
MIN_FREQ=0.01

CHROMOSOMES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/analyses/popmaps/chromosomes.txt
IDENTIFIER=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $CHROMOSOMES)

python /storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/UTILS/Pop_specific_SNP_selector.py /storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/analyses/cohort_1/freq/ ${IDENTIFIER}.frq $POPS $MIN_FREQ

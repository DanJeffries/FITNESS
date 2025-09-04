#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --array=1-5
#SBATCH --job-name=Freqs
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Analysis/vcftools/0.1.15

FILTERED_VCF_DIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/Filtered_VCFs/
N_MISSING_INFO=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/scripts/N_missing_info.txt
POP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $N_MISSING_INFO | cut -f2)
VCF_PATH=$FILTERED_VCF_DIR/McGill_filtered_masked_${POP}.vcf.gz.recode.vcf
TEMP_DIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/

vcftools --vcf $VCF_PATH \
         --out ${FILTERED_VCF_DIR}/${POP}_frequencies \
         --freq \
	 --temp $TEMP_DIR 




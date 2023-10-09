#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --array=1-5
#SBATCH --job-name=VCF_filter_0.6
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Analysis/vcftools/0.1.15

JOINT_CALL_VCF=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/Unfiltered_VCF/Joint_called.vcf
FILTERED_VCF_DIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/Filtered_VCFs/
N_MISSING_INFO=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/scripts/N_missing_info.txt

POP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $N_MISSING_INFO | cut -f2)
N_MISSING=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $N_MISSING_INFO | cut -f3)

## Pop numbers and number allowed missing

#      N POP N_MISS	
#      6 CN  2
#     25 FG  5
#     28 LG  5
#     27 SL  5
#     16 TN  3
#     29 WB  5
#     25 WK  5
#     26 WT  5


if [ ! -d '$(FILTERED_VCF_DIR)' ]; then
   mkdir $FILTERED_VCF_DIR
fi

ALL_SAMPLE_NAMES=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/scripts/sample_names.txt
POP_SAMPLE_MAP=${FILTERED_VCF_DIR}/sample_names_${POP}.txt
MASK=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/scripts/repeat_mask_sorted.bed


grep $POP $ALL_SAMPLE_NAMES > $POP_SAMPLE_MAP

TEMP_DIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/

vcftools --gzvcf $JOINT_CALL_VCF \
         --out ${FILTERED_VCF_DIR}/McGill_filtered_masked_mac_${POP}.vcf.gz \
         --recode \
         --recode-INFO-all \
	 --keep $POP_SAMPLE_MAP \
         --remove-indels \
         --max-alleles 2 \
	 --mac 1 \
         --max-missing-count $N_MISSING \
         --minQ 30 \
         --minGQ 30 \
	 --exclude-bed $MASK \
	 --temp $TEMP_DIR# 




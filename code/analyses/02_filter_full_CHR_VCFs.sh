#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --export=NONE
#SBATCH --array=1 #-22
#SBATCH --job-name=VCFtools_filter
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load Anaconda3
eval "$(conda shell.bash hook)"

WD=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS

VCF_DIR=$WD/DV_calling/CHR_VCFs_fill_missing

VCF=$VCF_DIR/chr_${SLURM_ARRAY_TASK_ID}.vcf.gz

FILTERED_VCF_OUT_PREFIX=$WD/DV_calling/CHR_VCFs_fill_missing_FILTERED/chr_${SLURM_ARRAY_TASK_ID}

vcftools --gzvcf $VCF \
	 --maf 0.1 \
	 --max-maf 0.4 \
	 --recode \
	 --out $FILTERED_VCF_OUT_PREFIX

gzip $FILTERED_VCF_OUT_PREFIX.*  ## the "." is important here


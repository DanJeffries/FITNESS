#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=48G
#SBATCH --export=NONE
#SBATCH --array=1-9
#SBATCH --job-name=VCFtools_filter
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load Anaconda3
eval "$(conda shell.bash hook)"

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS

POPULATIONS="FG  LG  SL  SR  TL  WB  WK  WT mixed"
POP=$(echo $POPULATIONS | cut -f $SLURM_ARRAY_TASK_ID -d' ')
POPMAP=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/analyses/cohort_1_pops/${POP}.pop

POP_SUBSAMPLE_GZVCF_OUT=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/analyses/cohort_1_tests/${POP}_cohort_1.subsample_0.01.vcf.gz

POP_SUBSAMPLE_GZVCF_VARONLY_OUT_PREFIX=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/analyses/cohort_1_tests/${POP}_cohort_1.subsample_0.01_filtered

vcftools --gzvcf $POP_SUBSAMPLE_GZVCF_OUT \
	 --minDP 5 \
	 --maxDP 50 \
	 --max-missing 0.75 \
	 --maf 0.05 \
	 --max-maf 0.95 \
	 --recode \
	 --out $POP_SUBSAMPLE_GZVCF_VARONLY_OUT_PREFIX

gzip $POP_SUBSAMPLE_GZVCF_VARONLY_OUT_PREFIX


#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --array=1
#SBATCH --job-name=subsample_var_only
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load Anaconda3
eval "$(conda shell.bash hook)"

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS

POPULATIONS="FG  LG  SL  SR  TL  WB  WK  WT mixed"
POP=$(echo $POPULATIONS | cut -f $SLURM_ARRAY_TASK_ID -d' ')
POPMAP=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/analyses/cohort_1_pops/${POP}.pop
POP_GZVCF=$WD/analyses/cohort_1_tests/${POP}_cohort_1.unfiltered_BCFtools.vcf.gz

## to make the test VCF I will downsample the FG population again
#POP_GZVCF=$WD/analyses/cohort_1_tests/FG_cohort_1.subsample_0.01.vcf.gz

r=0.01

POP_SUBSAMPLE_GZVCF_OUT=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/analyses/cohort_1_tests/${POP}_cohort_1.subsample_${r}.vcf.gz

#POP_SUBSAMPLE_GZVCF_OUT=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/analyses/cohort_1_tests/${POP}_cohort_1.subsample_0.0001.vcf.gz

conda activate vcflib_env

bcftools view $POP_GZVCF |\
\
vcfrandomsample -r $r |\
\
bcftools view -O z > $POP_SUBSAMPLE_GZVCF_OUT 





#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --array=9
#SBATCH --job-name=BCFtools_split
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load BCFtools/1.12-GCC-10.3.0

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS
GZVCF=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_calling/concatenated_cohort_VCFs/cohort_1_concatenated.vcf.gz

## First, split the VCF down into one for a specific population (doing this in an array)

POPULATIONS="FG  LG  SL  SR  TL  WB  WK  WT mixed"
POP=$(echo $POPULATIONS | cut -f $SLURM_ARRAY_TASK_ID -d' ')
POPMAP=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/analyses/cohort_1_pops/${POP}.pop

POP_VCF_OUT=$WD/analyses/cohort_1_tests/${POP}_cohort_1.unfiltered_BCFtools.vcf.gz

bcftools view $GZVCF \
	      -S $POPMAP \
	      --force-samples \
     	      -O z \
     	      > $POP_VCF_OUT



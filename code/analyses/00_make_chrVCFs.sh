#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --export=NONE
#SBATCH --array=1-22 # 1 job per chromosome
#SBATCH --job-name=Merge_chr_VCFs
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

WD=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/

GZVCFs=$(ls ${WD}/analyses/per_cohort/filtered/cohort_*/chr_${SLURM_ARRAY_TASK_ID}.filtered.recode.vcf.gz)

CHR_VCF_OUTDIR=$WD/DV_calling/CHR_VCFs_prefiltered

if [ ! -d "$CHR_VCF_OUTDIR" ]
then
    mkdir $CHR_VCF_OUTDIR
fi

CHR_VCF=$CHR_VCF_OUTDIR/chr_${SLURM_ARRAY_TASK_ID}.allsamples.prefiltered.vcf.gz

bcftools merge -O z \
	       -o $CHR_VCF \
	       $GZVCFs 



#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --export=NONE
#SBATCH --array=1-220
#SBATCH --job-name=VCFtools_filter
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load Anaconda3
eval "$(conda shell.bash hook)"

WD=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS

ITER_FILE=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/analyses/cohort_chrom_iterfile.txt

COHORT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE | cut -f1)
CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE | cut -f2)

VCF=$WD/DV_calling/GLnexus/cohort_${COHORT}/cohort_${COHORT}_chr_${CHROM}.vcf.gz

STAT=missing-indv ## used to make the output dir and to pass as an arg to vcftools

FILTERED_OUTDIR=$WD/analyses/per_cohort/filtered/cohort_${COHORT}

if [ ! -d "$FILTERED_OUTDIR" ]; then
        mkdir -p $FILTERED_OUTDIR
fi

FILTERED_VCF_OUT_PREFIX=$FILTERED_OUTDIR/chr_${CHROM}.filtered

#vcftools --gzvcf $VCF \
#	 --minDP 5 \
#	 --maxDP 75 \
#	 --mac 2 \
#	 --max-missing 0.25 \
#	 --recode \
#	 --out $FILTERED_VCF_OUT_PREFIX

#gunzip $FILTERED_VCF_OUT_PREFIX.*  ## the "." is important here

#bgzip $FILTERED_VCF_OUT_PREFIX.*

tabix $FILTERED_VCF_OUT_PREFIX.*

#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=01:00:00 
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --array=1-220
#SBATCH --job-name=VCFtools_STATS_ChrVCFsALL
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load VCFtools/0.1.16-GCC-10.3.0

WD=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS

ITER_FILE=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/analyses/cohort_chrom_iterfile.txt

COHORT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE | cut -f1)
CHROM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE | cut -f2)

VCF=$WD/DV_calling/GLnexus/cohort_${COHORT}/cohort_${COHORT}_chr_${CHROM}.vcf.gz

STAT=missing-indv ## used to make the output dir and to pass as an arg to vcftools

## potential stats to use are below, can swap them out of the command as needed (only one per command is allowed)

#        --het \
#        --hardy \
#        --freq2
#        --SNPdensity \
#	 --missing-indv \
#        --missing-site \
#        --site-mean-depth \

## make the output directory
STATS_DIR=$WD/analyses/per_cohort/missing_indv

if [ ! -d "$STATS_DIR" ]; then
	mkdir -p $STATS_DIR
fi

## calculate statistics

vcftools --gzvcf $VCF \
	 --${STAT} \
	 --out $STATS_DIR/cohort_${COHORT}_chr_${CHROM}


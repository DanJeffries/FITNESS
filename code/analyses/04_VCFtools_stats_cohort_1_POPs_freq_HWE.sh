#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=04:00:00 
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --export=NONE
#SBATCH --array=1-22
#SBATCH --job-name=VCFtools_cohort_1_FREQ
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load VCFtools/0.1.16-GCC-10.3.0

WD=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS

#FG  LG  SL  SR  TL  WB  WK  WT
POP=$1

VCF=$WD/DV_calling/GLnexus/cohort_1/cohort_1_chr_${SLURM_ARRAY_TASK_ID}.vcf.gz
POPMAP=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/analyses/popmaps/corhort_1/${POP}.pop

## potential stats to use are below, can swap them out of the command as needed (only one per command is allowed)

#        --het \
#        --hardy \
#        --freq2
#        --SNPdensity \
#	 --missing-indv \
#        --missing-site \
#        --site-mean-depth \
#	 --extract-FORMAT-info AD \

## make the output directory
STATS_OUTDIR=$WD/analyses/cohort_1/freq/$POP

if [ ! -d "$STATS_OUTDIR" ]; then
	mkdir -p $STATS_OUTDIR
fi

## calculate statistics

vcftools --gzvcf $VCF \
	 --keep $POPMAP \
	 --freq \
	 --remove-indels \
	 --max-missing 0.2 \
	 --out $STATS_OUTDIR/${POP}_chr_${SLURM_ARRAY_TASK_ID}


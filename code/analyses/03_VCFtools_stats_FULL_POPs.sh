#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=02:00:00 
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --export=NONE
#SBATCH --array=1-22
#SBATCH --job-name=VCFtools_STATS_FULL_PERPOP_PERCHROM
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load VCFtools/0.1.16-GCC-10.3.0

WD=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS

#FG  LG  SL  SR  TL  WB  WK  WT
POP="WT"

VCF=$WD/DV_calling/CHR_VCFs_prefiltered_NOFILL/chr_${SLURM_ARRAY_TASK_ID}.allsamples.prefiltered.vcf.gz
POPMAP=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/analyses/popmaps/full/${POP}_samples.txt

STAT="extract-FORMAT-info" ## used to make the output dir and to pass as an arg to vcftools

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
STATS_DIR=$WD/analyses/full_prefiltered_NOFILL_perpop/$STAT/$POP

if [ ! -d "$STATS_DIR" ]; then
	mkdir -p $STATS_DIR
fi

## calculate statistics

vcftools --gzvcf $VCF \
	 --${STAT} AD \
	 --keep $POPMAP \
	 --out $STATS_DIR/pop_${POP}_chr_${SLURM_ARRAY_TASK_ID}.prefiltered.filled


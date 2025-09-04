#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=01:00:00 
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --array=1-8
#SBATCH --job-name=VCFtools_POP_stats_females
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load VCFtools/0.1.16-GCC-10.3.0

WD=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/analyses/chr19_females_perpop

POPULATIONS="FG  LG  SL  SR  TL  WB  WK  WT"
POP=$(echo $POPULATIONS | cut -f $SLURM_ARRAY_TASK_ID -d' ')
GZVCF=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_calling/chr19_sepsexes/female/all_samples_chr_19_females.recode.vcf.gz
POPMAP=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/analyses/popmaps/sepsexes_chr19/females/${POP}_females.txt
_

if [ ! -d "$STATS_DIR" ]; then
	mkdir $STATS_DIR
fi

## calculate heterozygosities

vcftools --gzvcf $GZVCF \
	 --hardy \
	 --out $WD/${POP}_stats \
	 --keep $POPMAP \
#	 --het \
# 	 --freq2
#	 --SNPdensity \
#	 --missing-site \
#	 --site_depth 

## calculate SNPdensity



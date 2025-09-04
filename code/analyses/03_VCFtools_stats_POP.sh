#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=01:00:00 
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --array=1-8
#SBATCH --job-name=VCFtools_POP_stats
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load VCFtools/0.1.16-GCC-10.3.0

WD=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_calling/concatenated_cohort_VCFs/best_400

POPULATIONS="FG  LG  SL  SR  TL  WB  WK  WT"
POP=$(echo $POPULATIONS | cut -f $SLURM_ARRAY_TASK_ID -d' ')
GZVCF=$WD/cohort_1_best400.recode.vcf.gz
POPMAP=$WD/popmaps/${POP}_top50.txt
_

if [ ! -d "$STATS_DIR" ]; then
	mkdir $STATS_DIR
fi

## calculate heterozygosities

vcftools --gzvcf $GZVCF \
	 --freq2 \
	 --out $WD/${POP}_stats \
	 --keep $POPMAP \
#	 --het \
# 	 --freq2
#	 --SNPdensity \
#	 --missing-site \
#	 --site_depth 

## calculate SNPdensity



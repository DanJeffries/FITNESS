#!/bin/bash

#SBATCH --partition=epyc2
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

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS

POPULATIONS="FG  LG  SL  SR  TL  WB  WK  WT"
POP=$(echo $POPULATIONS | cut -f $SLURM_ARRAY_TASK_ID -d' ')
POP_GZVCF=$WD/analyses/cohort_1_tests/${POP}_cohort_1.subsample_0.01.vcf.gz

STATS_DIR=$WD/analyses/cohort_1_tests/stats

if [ ! -d "$STATS_DIR" ]; then
	mkdir $STATS_DIR
fi

## calculate heterozygosities

vcftools --gzvcf $POP_GZVCF \
	 --het \
	 --out $STATS_DIR/${POP}_stats
# 	  --freq2
#	 --SNPdensity \
#	 --missing-site \
#	 --site_depth 

## calculate SNPdensity



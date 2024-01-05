#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --export=NONE
#SBATCH --array=1-92
#SBATCH --job-name=SAMSTATS
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Analysis/samtools/1.10

ID=$SLURM_ARRAY_TASK_ID
ITER_FILE=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/bamMDRG/bams.txt
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE | cut -f1 -d'.')

BAM_DIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/bamMDRG/
BAM_IN=${BAM_DIR}/${SAMPLE_NAME}.MD.RG.bam

STATS_DIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/bamMDRG/SamStats/
STATS_OUT=${STATS_DIR}/${SAMPLE_NAME}.stats

if [ ! -d "$STATS_DIR" ]; then
   mkdir $STATS_DIR
fi

samtools stats $BAM_IN > $STATS_OUT



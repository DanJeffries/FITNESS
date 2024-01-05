#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --export=NONE
#SBATCH --array=1-92
#SBATCH --job-name=SAMCOV
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Analysis/samtools/1.10

ID=$SLURM_ARRAY_TASK_ID
ITER_FILE=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/bamMDRG/bams.txt
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE | cut -f1 -d'.')

BAM_DIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/bamMDRG/
BAM_IN=${BAM_DIR}/${SAMPLE_NAME}.MD.RG.bam

DEPTHS_DIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/SamCov/
DEPTHS_OUT=${DEPTHS_DIR}/${SAMPLE_NAME}.cov

if [ ! -d "$DEPTHS_DIR" ]; then
   mkdir $DEPTHS_DIR
fi

samtools coverage $BAM_IN > $DEPTHS_OUT

## Calculate some averages and add to the end of each file 
egrep -v '#|chrUn|chrM' $DEPTHS_OUT | awk '{ sum += $6; n++ } END { if (n > 0) print "\n## Avg. perc. bases. covered = " sum / n; }' >> $DEPTHS_OUT
egrep -v '#|chrUn|chrM' $DEPTHS_OUT | awk '{ sum += $7; n++ } END { if (n > 0) print "## Avg. coverage = " sum / n; }' >> $DEPTHS_OUT


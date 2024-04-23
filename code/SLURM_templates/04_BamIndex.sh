#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60G
#SBATCH --export=NONE
#SBATCH --array=1-512
#SBATCH --job-name=BamIndex
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Aligner/bwa/0.7.17
module load UHTS/Analysis/samblaster/0.1.24 
module load UHTS/Analysis/samtools/1.10
module load UHTS/Analysis/sambamba/0.7.1

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS_DATA/
DATA_DIR=$WD/trimmed
ITER_FILE=$DATA_DIR/sample_names.txt
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE)

OUTDIR=$WD/bams/

BAM_FIX_COORDSORT=${OUTDIR}/${SAMPLE_NAME}.fixmate.coordsorted.bam

samtools index -@ 20 $BAM_FIX_COORDSORT

#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --export=NONE
#SBATCH --array=1-18
#SBATCH --job-name=FASTQC
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Quality_control/fastqc/0.11.9

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training

ID=$SLURM_ARRAY_TASK_ID
ITERFILE=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/code/samples.txt

DATA_DIR=${WD}/trimmed
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITERFILE)

OUTDIR=${WD}/trimmed/fastqc

if [ ! -d "$OUTDIR" ]
then
    mkdir -p $OUTDIR
fi

fastqc $WD/trimmed/${SAMPLE_NAME}.R1.trimmed.fastq.gz -t 6 -o $OUTDIR
fastqc $WD/trimmed/${SAMPLE_NAME}.R2.trimmed.fastq.gz -t 6 -o $OUTDIR

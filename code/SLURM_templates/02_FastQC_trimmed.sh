#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --export=NONE
#SBATCH --array=21-1020
#SBATCH --job-name=FASTQC
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Quality_control/fastqc/0.11.9

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS

ID=$SLURM_ARRAY_TASK_ID
ITERFILE=${WD}/data/trimmed/trimmed_fastqs.txt

DATA_DIR=${WD}/data/trimmed
FASTQ=${DATA_DIR}/$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITERFILE)

OUTDIR=${WD}/data/FastQC_trimmed

if [ ! -d "$OUTDIR" ]
then
    mkdir -p $OUTDIR
fi

fastqc $FASTQ -t 6 -o $OUTDIR

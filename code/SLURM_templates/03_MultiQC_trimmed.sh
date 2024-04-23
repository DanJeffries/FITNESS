#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --export=NONE
#SBATCH --job-name=MultiQC_trim
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load UHTS/Analysis/MultiQC/1.8

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS
DATA_DIR=${WD}/data/FastQC_trimmed
OUT_DIR=${WD}/data/MultiQC_trimmed

if [ ! -d "$OUT_DIR" ]
then
    mkdir -p $OUT_DIR
fi  

multiqc $DATA_DIR -o $OUT_DIR


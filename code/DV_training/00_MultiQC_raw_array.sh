#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --export=NONE
#SBATCH --job-name=MultiQC
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load UHTS/Analysis/MultiQC/1.8

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/

RUN_DATA_DIR=$WD/raw/fastqc
RUN_OUT_DIR=$WD/MultiQC_raw

if [ ! -d "$RUN_OUT_DIR" ]
then
    mkdir -p $RUN_OUT_DIR
fi  

multiqc --ignore *_I[1-2]_* $RUN_DATA_DIR -o $RUN_OUT_DIR


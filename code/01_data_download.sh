#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --export=NONE
#SBATCH --array=1-100%6
#SBATCH --job-name=FITNESS_DOWNLOAD
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load UHTS/Analysis/sratoolkit/2.10.7

ID=$SLURM_ARRAY_TASK_ID
URLS=/storage/scratch/iee/dj20y461/FITNESS/TEST_12_22/raw/URLs.txt
URL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $URLS)
DATA_DIR=/storage/scratch/iee/dj20y461/FITNESS/TEST_12_22/raw/

cd $DATA_DIR

wget $URL

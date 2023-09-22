#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --export=NONE
#SBATCH --array=1-100%6

60,63,64,127,128,129,130,131,132,133,134,135,136,137,138,139,140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,161,162,163,164,165,166,167,168,169,170,171,172,173,174,175,176,177,178,179%6
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

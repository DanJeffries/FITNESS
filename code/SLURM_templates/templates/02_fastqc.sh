#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --export=NONE
#SBATCH --array=139,156,161,318,331,336,352,355
#SBATCH --job-name=FASTQC
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Quality_control/fastqc/0.11.9

ID=$SLURM_ARRAY_TASK_ID
DATA_DIR=/storage/scratch/iee/dj20y461/FITNESS/TEST_12_22/raw/

ls ${DATA_DIR}*gz > $DATA_DIR/fastq_files.txt

FASTQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $DATA_DIR/fastq_files.txt)

fastqc $FASTQ -t 6

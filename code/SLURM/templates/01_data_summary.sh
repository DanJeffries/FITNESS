#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --export=NONE
#SBATCH --array=9
#SBATCH --job-name=R2_DATA_SUM
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load UHTS/Analysis/sratoolkit/2.10.7

RUN=2
ID=$SLURM_ARRAY_TASK_ID
FASTQs=/storage/scratch/iee/dj20y461/FITNESS/data/raw/fastq_links/Run${RUN}_R1_paths.txt
FASTQ=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $FASTQs)
OUTFILE=/storage/scratch/iee/dj20y461/FITNESS/plot_data/read_counts_$RUN.tsv

READS=$(expr $(zcat $FASTQ | wc -l) / 4)

printf "%s\t%s\n" "$FASTQ" "$READS" >> $OUTFILE


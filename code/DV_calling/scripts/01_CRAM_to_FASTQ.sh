#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=04:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=20G
#SBATCH --export=NONE
#SBATCH --array=2-340%8
#SBATCH --job-name=CRAM2FASTQ
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

CRAM_FILES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_calling/scripts/cram_files.txt
CRAM_DIR=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/raw_BACKUP/FITNESS_Run_S/

CRAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $CRAM_FILES)

CRAM_NAME=$(echo $CRAM | rev | cut -f1 -d'/' | rev)

FASTQ_NAMES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_calling/scripts/CRAM_SAMPLE_NAMES.txt

FASTQ_NAME=$(grep "$CRAM_NAME" $FASTQ_NAMES | cut -f2)

echo "converting $CRAM_NAME to $FASTQ_NAME"

samtools view -b $CRAM -@ 20| samtools sort -n | samtools fastq -@ 20 -1 $CRAM_DIR/${FASTQ_NAME}_2-0000000_S00_L000_R1_000.fastq -2 $CRAM_DIR/${FASTQ_NAME}_2-0000000_S00_L000_R2_000.fastq -

gzip $CRAM_DIR/${FASTQ_NAME}_2-0000000_S00_L000_R1_000.fastq
gzip $CRAM_DIR/${FASTQ_NAME}_2-0000000_S00_L000_R2_000.fastq


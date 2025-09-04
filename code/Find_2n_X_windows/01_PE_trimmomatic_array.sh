#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --export=NONE
#SBATCH --array=1-25
#SBATCH --job-name=TRIM
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Analysis/trimmomatic/0.36 

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/Find_2n_X_windows

ID=$SLURM_ARRAY_TASK_ID
ITER_FILE=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/Find_2n_X_windows/samples_new.txt
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE)

ADAPTERS=/software/UHTS/Analysis/trimmomatic/0.36/adapters/NexteraPE-PE.fa

RAW_DIR=$WD/raw_links/new/
TRIMMED_OUTDIR=$WD/trimmed

if [ ! -d "$TRIMMED_OUTDIR" ] 
then
    mkdir $TRIMMED_OUTDIR
fi

## OUT FILES

SAMPLE_1_kept=$(echo ${SAMPLE_NAME}.R1.trimmed.fastq.gz)
SAMPLE_1_unpaired=$(echo ${SAMPLE_NAME}.R1.unpaired.fastq.gz)
SAMPLE_2_kept=$(echo ${SAMPLE_NAME}.R2.trimmed.fastq.gz)
SAMPLE_2_unpaired=$(echo ${SAMPLE_NAME}.R2.unpaired.fastq.gz)

trimmomatic PE -threads 4 \
	$RAW_DIR/$SAMPLE_NAME*_R1_*fastq.gz $RAW_DIR/$SAMPLE_NAME*_R2_*fastq.gz \
	${TRIMMED_OUTDIR}/${SAMPLE_1_kept} ${TRIMMED_OUTDIR}/${SAMPLE_1_unpaired} \
	${TRIMMED_OUTDIR}/${SAMPLE_2_kept} ${TRIMMED_OUTDIR}/${SAMPLE_2_unpaired} \
	ILLUMINACLIP:${ADAPTERS}:3:30:10 SLIDINGWINDOW:4:20 MINLEN:25 HEADCROP:7


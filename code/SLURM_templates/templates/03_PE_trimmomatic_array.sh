#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --export=NONE
#SBATCH --array=1-92
#SBATCH --job-name=TRIM
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Analysis/trimmomatic/0.36 

ID=$SLURM_ARRAY_TASK_ID
ITER_FILE=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/scripts/sample_names.txt
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE)
DATA_DIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/raw_concat/
ADAPTERS=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/scripts/illumina_adapters.fa

TRIMMED_OUTDIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/trimmed/

if [ ! -d "$TRIMMED_OUTDIR" ] 
then
    mkdir $TRIMMED_OUTDIR
fi

## OUT FILES

SAMPLE_1_kept=$(echo ${SAMPLE_NAME}_1.trimmed.fastq.gz)
SAMPLE_1_unpaired=$(echo ${SAMPLE_NAME}_1.unpaired.fastq.gz)
SAMPLE_2_kept=$(echo ${SAMPLE_NAME}_2.trimmed.fastq.gz)
SAMPLE_2_unpaired=$(echo ${SAMPLE_NAME}_2.unpaired.fastq.gz)

trimmomatic PE -threads 4 \
	${DATA_DIR}/${SAMPLE_NAME}.R1.fastq.gz ${DATA_DIR}/${SAMPLE_NAME}.R2.fastq.gz \
	${TRIMMED_OUTDIR}/${SAMPLE_1_kept} ${TRIMMED_OUTDIR}/${SAMPLE_1_unpaired} \
	${TRIMMED_OUTDIR}/${SAMPLE_2_kept} ${TRIMMED_OUTDIR}/${SAMPLE_2_unpaired} \
	ILLUMINACLIP:${ADAPTERS}:3:30:10 SLIDINGWINDOW:4:20 MINLEN:25 HEADCROP:7


#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=03:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --export=NONE
#SBATCH --array=511-765 ## Run2
#SBATCH --job-name=TRIM
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Analysis/trimmomatic/0.36 

ID=$SLURM_ARRAY_TASK_ID
ITER_FILE=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/data/raw/Unique_R1_raw_paths.txt

R1_FASTQ_PATH=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE)
R2_FASTQ_PATH=$(echo $R1_FASTQ_PATH | sed 's/_R1_/_R2_/g')
SAMPLE_NAME=$(echo $R1_FASTQ_PATH | rev | cut -d '/' -f 1  | rev | cut -f-4 -d'_')

ADAPTERS=/software/UHTS/Analysis/trimmomatic/0.36/adapters/NexteraPE-PE.fa

TRIMMED_OUTDIR=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/data/trimmed

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
	$R1_FASTQ_PATH $R2_FASTQ_PATH \
	${TRIMMED_OUTDIR}/${SAMPLE_1_kept} ${TRIMMED_OUTDIR}/${SAMPLE_1_unpaired} \
	${TRIMMED_OUTDIR}/${SAMPLE_2_kept} ${TRIMMED_OUTDIR}/${SAMPLE_2_unpaired} \
	ILLUMINACLIP:${ADAPTERS}:3:30:10 SLIDINGWINDOW:4:20 MINLEN:25 HEADCROP:7


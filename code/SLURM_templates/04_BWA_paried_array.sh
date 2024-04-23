#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --export=NONE
#SBATCH --array=1-92
#SBATCH --job-name=BWA
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Aligner/bwa/0.7.17

ID=$SLURM_ARRAY_TASK_ID
ITER_FILE=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/scripts/sample_names.txt
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE)
TRIMMED_DIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/trimmed/
SAM_DIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/sam/

if [ ! -d "$SAM_DIR" ]; then
   mkdir $SAM_DIR
fi

GENOME_IDX=/storage/scratch/iee/dj20y461/Gacu_assembly/No_Y/BWA/stickleback_v5_assembly_NoY.fa

SAMPLE_1=$(echo ${TRIMMED_DIR}/${SAMPLE_NAME}_1.trimmed.fastq.gz)
SAMPLE_2=$(echo ${TRIMMED_DIR}/${SAMPLE_NAME}_2.trimmed.fastq.gz)

SAM=${SAM_DIR}/${SAMPLE_NAME}.sam

bwa mem -t 8 -o $SAM $GENOME_IDX $SAMPLE_1 $SAMPLE_2

echo "Reads aligned"

### Sort, convert to bam and index

BAM_DIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/bam/

if [ ! -d "$BAM_DIR" ]; then
   mkdir $BAM_DIR
fi

<<<<<<< HEAD
BAM_FIX_COORDSORT=${OUTDIR}/${SAMPLE_NAME}.fixmate.coordsorted.bam


#bwa mem -t 20 \
#        -R $(echo "@RG\tID:${SAMPLE_NAME}\tPL:ILLUMINA\tPU:${FLOWCELL}.${LANE}\tSM:${SAMPLE_NAME}\tCN:MCGILL\tBC:${BARCODES}") \
#        $GENOME_IDX \
#        $SAMPLE_R1 \
#        $SAMPLE_R2 | \
#samblaster | \
#samtools fixmate \
#	-m \
#        -@ 20 \
#        - \
#	- | \
#samtools sort \
#	-@ 20 \
#	-O BAM \
#	-o $BAM_FIX_COORDSORT \
#	-

samtools index -c -@ 20 $BAM_FIX_COORDSORT
=======
BAM=${BAM_DIR}/${SAMPLE_NAME}.bam

# SORT AND OUPUT BAM

samtools sort -@ 16 -o $BAM $SAM

echo "BAM Sorted" 

rm $SAM

echo "SAM removed"

## INDEX

samtools index -c -@ 16 $BAM

echo "BAM indexed"
>>>>>>> parent of 6feef65... BWA script

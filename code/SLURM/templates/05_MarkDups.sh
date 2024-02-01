#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --export=NONE
#SBATCH --array=57,89
#SBATCH --job-name=MDRG
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Analysis/samtools/1.10
module load UHTS/Analysis/picard-tools/2.21.8


ID=$SLURM_ARRAY_TASK_ID
ITER_FILE=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/scripts/sample_names.txt
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE)

BAM_DIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/bam/
BAM_IN=${BAM_DIR}/${SAMPLE_NAME}.bam

# Mark Duplicates

MARKED_DUPS_OUTDIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/bamMD/
MARKED_DUPS_OUT_BAM=${MARKED_DUPS_OUTDIR}/${SAMPLE_NAME}.MD.bam
MARKED_DUPS_OUT_METRICS=${MARKED_DUPS_OUTDIR}/${SAMPLE_NAME}.MDmetrics.txt

if [ ! -d "$MARKED_DUPS_OUTDIR" ]; then
   mkdir $MARKED_DUPS_OUTDIR
fi

TMP_DIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/javatmp/

if [ ! -d "$TMP_DIR" ]; then
   mkdir $TMP_DIR
fi

picard-tools MarkDuplicates \
      I=${BAM_IN} \
      O=${MARKED_DUPS_OUT_BAM} \
      M=${MARKED_DUPS_OUT_METRICS} \
      OPTICAL_DUPLICATE_PIXEL_DISTANCE=2500 \
      REMOVE_DUPLICATES=true \
      TMP_DIR=$TMP_DIR

MD_status=$?

if [ $MD_status -ne 0 ]; then
    echo "picard-tools MarkDuplicates failed"
    exit $MD_status
fi


echo "MarkDuplicates passed" 

# Remove original bam

#rm $BAM_IN

# Index new bam

samtools index -c ${MARKED_DUPS_OUT_BAM}

echo "Dups marked and bam indexed . . ."

# Add Readgroups

#RG_OUTDIR=/storage/scratch/iee/dj20y461/PopHistory/Global/GE/bamMDRG/
#RG_OUT_BAM=${RG_OUTDIR}/${SAMPLE_NAME}.MD.RG.bam
#
#if [ ! -d "$RG_OUTDIR" ]; then
#   mkdir $RG_OUTDIR
#fi
#
#picard-tools AddOrReplaceReadGroups \
#       I=$MARKED_DUPS_OUT_BAM  \
#       O=${RG_OUT_BAM} \
#       RGID=GE \
#       RGLB=lib_GE \
#       RGPL=illumina \
#       RGPU=unit_GE \
#       RGSM=${SAMPLE_NAME}
#
#samtools index -c ${RG_OUT_BAM}
#
#echo "New bam indexed . . . All done!"


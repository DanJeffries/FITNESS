#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=80G
#SBATCH --export=NONE
#SBATCH --array=1-92
#SBATCH --job-name=AddRG
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Analysis/samtools/1.10
module load UHTS/Analysis/picard-tools/2.21.8


ID=$SLURM_ARRAY_TASK_ID
ITER_FILE=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/scripts/sample_names.txt
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE)

MARKED_DUPS_OUTDIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/bamMD/
MARKED_DUPS_OUT_BAM=${MARKED_DUPS_OUTDIR}/${SAMPLE_NAME}.MD.bam

TMP_DIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/javatmp/


# Add Readgroups

RG_OUTDIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/bamMDRG/
RG_OUT_BAM=${RG_OUTDIR}/${SAMPLE_NAME}.MD.RG.bam

if [ ! -d "$RG_OUTDIR" ]; then
   mkdir $RG_OUTDIR
fi

picard-tools AddOrReplaceReadGroups \
       I=$MARKED_DUPS_OUT_BAM  \
       O=${RG_OUT_BAM} \
       RGID=RG1 \
       RGLB=lib_1 \
       RGPL=illumina \
       RGPU=unit_1 \
       RGSM=${SAMPLE_NAME} \
       TMP_DIR=$TMP_DIR

samtools index -c ${RG_OUT_BAM}

echo "New bam indexed . . . All done!"


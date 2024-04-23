#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=48G
#SBATCH --export=NONE
#SBATCH --array=2-187
#SBATCH --job-name=AddRG_MCGILL
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Analysis/samtools/1.10
module load UHTS/Analysis/picard-tools/2.21.8


ID=$SLURM_ARRAY_TASK_ID

# Add readgroup info

RG_INFO=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_01_23/scripts/RG_info.txt

SAMPLE_RG=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $RG_INFO | cut -f1)
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $RG_INFO | cut -f2)

MARKED_DUPS_OUTDIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_01_23/bamMD/
MARKED_DUPS_OUT_BAM=${MARKED_DUPS_OUTDIR}/${SAMPLE_NAME}.MD.bam

RG_OUTDIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_01_23/bamMDRG/
RG_OUT_BAM=${RG_OUTDIR}/${SAMPLE_NAME}.MDRG.bam

if [ ! -d "$RG_OUTDIR" ]; then
   mkdir $RG_OUTDIR
fi

picard-tools AddOrReplaceReadGroups \
       I=$MARKED_DUPS_OUT_BAM  \
       O=${RG_OUT_BAM} \
       RGID=$SAMPLE_RG \
       RGLB=MCGILL \
       RGPL=illumina \
       RGPU=$SAMPLE_RG \
       RGSM=${SAMPLE_NAME}

RG_status=$?

if [ $RG_status -ne 0 ]; then
    echo "picard-tools AddOrReplaceReadGroups failed"
    exit $RG_status
fi

echo "AddOrReplaceReadGroups passed"

samtools index -c ${RG_OUT_BAM}

echo "New bam indexed . . . All done!"

rm $MARKED_DUPS_OUT_BAM
rm ${MARKED_DUPS_OUT_BAM}.csi

echo "Old bam removed""


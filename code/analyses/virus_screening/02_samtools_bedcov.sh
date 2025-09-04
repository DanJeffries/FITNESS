#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --array=326-1012
#SBATCH --job-name=bedcov_virus_screening
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module add Trimmomatic/0.39-Java-11
module add BWA/0.7.17-GCC-10.3.0
module add samblaster/0.1.26-GCC-10.3.0
module add SAMtools/1.13-GCC-10.3.0
module add Sambamba/0.8.2-GCC-10.3.0

echo "start time"
date

WD=/storage/research/iee_temp_dj20y461/DV_calling/virus_screening

ID=$SLURM_ARRAY_TASK_ID
ITER_FILE=/storage/research/iee_temp_dj20y461/DV_calling/virus_screening/bams/bam_names.list

SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE)

BAMDIR=$WD/bams/

BAM_FIX_COORDSORT=${BAMDIR}/${SAMPLE_NAME}.fixmate.coordsorted.bam

IRIDI_BED_FILE=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/ref/Iridovirus/PQ335173_4_iridovirus_CDS.bed
HERPES_BED_FILE=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/ref/Herpesvirus/GasAcuHV-C3_sans_ORF7_8.bed
ECHO_BED_FILE=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/ref/Herpesvirus/GasAcuHV-C3_ORF7_8.bed

IRIDI_COV_OUT=${BAMDIR}/${SAMPLE_NAME}.iridovirus.bedcov
HERPES_COV_OUT=${BAMDIR}/${SAMPLE_NAME}.herpes_sans_ORF7_8.bedcov
ECHO_COV_OUT=${BAMDIR}/${SAMPLE_NAME}.herpesORF7_8.bedcov

echo $IRIDI_BED_FILE
samtools bedcov $IRIDI_BED_FILE $BAM_FIX_COORDSORT > $IRIDI_COV_OUT

echo $HERPES_BED_FILE
samtools bedcov $HERPES_BED_FILE $BAM_FIX_COORDSORT > $HERPES_COV_OUT

echo $ECHO_BED_FILE
samtools bedcov $ECHO_BED_FILE $BAM_FIX_COORDSORT > $ECHO_COV_OUT



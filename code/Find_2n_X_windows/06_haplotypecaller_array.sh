#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=10:00:00
#eBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --array=11-2100
#SBATCH --job-name=HapCall
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Analysis/GenomeAnalysisTK/4.1.3.0
module load UHTS/Analysis/picard-tools/2.9.0

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/Find_2n_X_windows/

ITERFILE=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/Find_2n_X_windows/sample_interval.map

SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITERFILE | cut -f1 -d' ')
INTERVAL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITERFILE | cut -f2 -d' ')

GENOME=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna
DICT=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.dict

if [ ! -f "$DICT" ]; then
   picard-tools CreateSequenceDictionary R=$GENOME O=$DICT
fi

BAM_DIR=$WD/bams
GVCF_DIR=$WD/GVCF

if [ ! -d "$GVCF_DIR" ]; then
   mkdir $GVCF_DIR
fi

BAM=${BAM_DIR}/${SAMPLE_NAME}.fixmate.coordsorted.bam

GVCF=${GVCF_DIR}/${SAMPLE_NAME}_${INTERVAL}.g.vcf.gz

java -jar $GATK_PATH/bin/GenomeAnalysisTK.jar HaplotypeCaller \
	-R $GENOME \
	-I $BAM \
	-L $INTERVAL \
	-ERC GVCF \
	--output $GVCF



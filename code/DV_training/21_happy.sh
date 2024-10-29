#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=02:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --export=NONE
#SBATCH --job-name=HAPPY_EVAL
#SBATCH --array=1 #-5
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

####script to run make_examples
WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/
REF=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_formatted_shortnames.fna
SAMPLES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/offspring.txt
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $SAMPLES)
CROSS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $SAMPLES | cut -f1 -d'_')

export APPTAINER_TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/apptained_tmp" #Set Singularity temporary dir
export TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/parralel_tmp" #Set global temporary dir for parallel

HAPPY_PATH=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV/happy.sif

# make output dir

if [ ! -d "$WD/evaluate/${SAMPLE}" ]; then
   mkdir -p $WD/evaluate/${SAMPLE}
fi

apptainer run \
--nv \
-B $WD:/home \
$HAPPY_PATH \
/opt/hap.py/bin/hap.py \
	/home/Filtered_variants/${SAMPLE}.ALL_TRUTH_VARS.vcf.gz \
	/home/test/${SAMPLE}_test_set.vcf.gz \
	-f /home/Confident_regions/${CROSS}_male_conf_regions_inc_vars.bed \
	-r /home/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_formatted_shortnames.fna \
	-o /home/evaluate/${SAMPLE}_test_happy_eval.output \
	--engine=vcfeval \
	--pass-only

#-l /home/training_regions/${CROSS}_test_partitions.bed


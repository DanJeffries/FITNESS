#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=02:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --export=NONE
#SBATCH --job-name=HAPPY_EVAL_GATK
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

####script to run make_examples
WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/
REF=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_formatted_shortnames.fna

export APPTAINER_TMPDIR="${WD}/apptainer_tmp" #Set Singularity temporary dir

if [ ! -d "$APPTAINER_TMPDIR" ]; then
   mkdir -p "$APPTAINER_TMPDIR"
fi

export TMPDIR="${WD}/parralel_tmp"

if [ ! -d "$TMPDIR" ]; then
   mkdir -p "$TMPDIR"
fi

HAPPY_PATH=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV/happy.sif

SAMPLE=FG_male_1
CROSS=FG

## Model to test
### this text file gives path to best checkpoint for each of the training runs being evaluated. Paths start from the path mounted as home in apptainer.
BEST_MODELS=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/20_testcalls_scripts/Step_1/best_checkpoints.txt

MODEL_DIR=GRIDSEARCH/gatk_baseline
CHECKPOINT=GATK

TEST_CALLSET=${MODEL_DIR}/FG_male_1.GATK.filtered_baseline.vcf.gz

## make outdir

HAPPY_OUTDIR=${MODEL_DIR}/happy_eval

if [ ! -d "$WD/$HAPPY_OUTDIR" ]; then
   mkdir -p "$WD/$HAPPY_OUTDIR"
fi

## all paths relative to $WD, which is the mount point in the container

TRUTH_SET=TRAINING_DATA/${SAMPLE}.CONF_VARS_ALL.vcf.gz
CONF_REGIONS=TRAINING_DATA/${SAMPLE}.CONF_REGIONS_MASKED.bed

TEST_PARTITIONS=$(awk '{print $1}' ${WD}/training_regions/${CROSS}_test_partitions.bed | paste -s -d , -)
echo $TEST_PARTITIONS
REF=ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_formatted_shortnames.fna

apptainer run \
--nv \
-B $WD:/home \
$HAPPY_PATH \
/opt/hap.py/bin/hap.py \
	/home/$TRUTH_SET \
	/home/$TEST_CALLSET \
	-f /home/$CONF_REGIONS \
	-r /home/$REF \
	-o /home/$HAPPY_OUTDIR/${SAMPLE}_${CHECKPOINT}_happy_eval.output \
	--engine=vcfeval \
	--pass-only \
	-l $TEST_PARTITIONS


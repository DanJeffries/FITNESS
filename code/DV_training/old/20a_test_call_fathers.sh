#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=06:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
##SBATCH --gres=gpu:rtx4090:1
#SBATCH --mem-per-cpu=10G
#SBATCH --export=NONE
#SBATCH --job-name=TRAIN_test_callFathers_CPU
#SBATCH --array=5
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

####script to run make_examples
WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/
REF=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_formatted_shortnames.fna
SAMPLES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/offspring.txt
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $SAMPLES)
CROSS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $SAMPLES | cut -f1 -d'_')
FATHER=${CROSS}_male_par

DV_PATH=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV/deepvariant_1.6.1-gpu.modified.sif
OPENBLAS_NUM_THREADS=1 #Set number of threads that OPENBLAS uses to avoid thread overflow error in numpy

## Model to test
MODEL=LR0.001_BS1024
CKPT=ckpt-1050
MODEL_SUBDIR=training_outs/TRAIN_ROUND1/$MODEL/checkpoints/$CKPT

# make output dir

OUTDIR=model_tests/$MODEL/calls/

if [ ! -d "$OUTDIR" ]; then
   mkdir -p "$OUTDIR"
fi

## set temp dirs for the run
export APPTAINER_TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/apptained_tmp/model_tests_$MODEL_$SAMPLE" #Set Singularity temporary dir
export TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/parralel_tmp/model_tests_$MODEL_$SAMPLE" #Set global temporary dir for parallel

if [ ! -d "$APPTAINER_TMPDIR" ]; then
   mkdir -p "$APPTAINER_TMPDIR"
fi

if [ ! -d "$TMPDIR" ]; then
   mkdir -p "$TMPDIR"
fi

apptainer run \
-B $WD:/home \
$DV_PATH \
/opt/deepvariant/bin/run_deepvariant \
--model_type WGS \
--customized_model "/home/${MODEL_SUBDIR}" \
--ref "/home/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_formatted_shortnames.fna" \
--reads "/home/bams/${FATHER}.fixmate.coordsorted.bam" \
--regions "/home/training_regions/${CROSS}_test_partitions.bed" \
--output_vcf "/home/test/${FATHER}_test_set.vcf.gz" \
--num_shards=20

#--nv \


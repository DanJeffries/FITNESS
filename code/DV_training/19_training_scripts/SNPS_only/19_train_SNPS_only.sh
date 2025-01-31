#!/bin/bash

#SBATCH --partition=gpu
#SBATCH --time=10:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:h100:1
#SBATCH --mem-per-cpu=50G
#SBATCH --export=NONE
#SBATCH --job-name=TRAIN_SNPS_only_test
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

####script to run make_examples
WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/
REF=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_formatted_shortnames.fna
SAMPLES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/offspring.txt
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $SAMPLES)

export APPTAINER_TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/apptained_tmp" #Set Singularity temporary dir
export TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/parralel_tmp" #Set global temporary dir for parallel

DV_PATH=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV/deepvariant_1.6.1-gpu.modified.sif
OPENBLAS_NUM_THREADS=1 #Set number of threads that OPENBLAS uses to avoid thread overflow error in numpy

# get step instructions

BS=64
LR=0.001
TUNE_EVERY=400

# make outdir 

OUTDIR=SNPS_only_training/

if [ ! -d "$WD/${OUTDIR}" ]; then
   mkdir -p $WD/${OUTDIR}/
fi


## Initial model to start training from


apptainer run \
--nv \
-B $WD:/home \
$DV_PATH \
/opt/deepvariant/bin/train \
--config=/home/dv_config.py:base \
--config.train_dataset_pbtxt="/home/examples_shuffled/train_SNPS_ONLY/examples_shuffled_config.pbtxt" \
--config.tune_dataset_pbtxt="/home/examples_shuffled/tune_SNPS_ONLY/tune_examples.SNPS_ONLY.config.pbtxt" \
--config.num_epochs=1 \
--config.learning_rate=$LR \
--config.num_validation_examples=0 \
--config.tune_every_steps=$TUNE_EVERY \
--experiment_dir=/home/${OUTDIR} \
--strategy=mirrored \
--config.batch_size=$BS



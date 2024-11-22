#!/bin/bash

#SBATCH --partition=gpu
#SBATCH --time=24:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:rtx4090:1
#SBATCH --mem-per-cpu=50G
#SBATCH --export=NONE
#SBATCH --job-name=TRAIN_lr_dropdown_test_0_01_to_0_0001
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

# make output dir

LR=0.0001
BS=1024
TUNE_STEPS=350

OUTDIR=training_outs/param_tests/manual_LR_dropdown/LR_0_01_to_LR${LR}_BS${BS}

if [ ! -d "$WD/${OUTDIR}" ]; then
   mkdir -p $WD/${OUTDIR}/
fi

apptainer run \
--nv \
-B $WD:/home \
$DV_PATH \
/opt/deepvariant/bin/train \
--config=/home/dv_config.py:base \
--config.train_dataset_pbtxt="/home/examples_shuffled/train/shuf3/examples_shuf3_sub_1_10_config.pbtxt" \
--config.tune_dataset_pbtxt="/home/examples_shuffled/tune/All_samples_tune_examples.dataset_config.pbtxt" \
--config.num_epochs=1 \
--config.learning_rate=$LR \
--config.num_validation_examples=0 \
--config.tune_every_steps=$TUNE_STEPS \
--experiment_dir=/home/${OUTDIR} \
--strategy=mirrored \
--config.batch_size=$BS \
--config.init_checkpoint="/home/training_outs/param_tests/test_0.01_1024/checkpoints/ckpt-500"




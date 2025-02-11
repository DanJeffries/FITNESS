#!/bin/bash

#SBATCH --partition=gpu
#SBATCH --time=24:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
##SBATCH --gres=gpu:h100:1
#SBATCH --gres=gpu:rtx4090:1
#SBATCH --mem-per-cpu=50G
#SBATCH --export=NONE
#SBATCH --array=5
#SBATCH --job-name=TRAIN_ALT	
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

####script to run make_examples
WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/

export APPTAINER_TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/apptained_tmp" #Set Singularity temporary dir
export TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/parralel_tmp" #Set global temporary dir for parallel

DV_PATH=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV/deepvariant_1.6.1-gpu.modified.sif
OPENBLAS_NUM_THREADS=1 #Set number of threads that OPENBLAS uses to avoid thread overflow error in numpy

# get step instructions

STEP_PARAMS="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/19_training_scripts/Step_1/step_parameters.txt"

STEP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $STEP_PARAMS | cut -f1)
RUN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $STEP_PARAMS | cut -f2)
BS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $STEP_PARAMS | cut -f3)
LR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $STEP_PARAMS | cut -f4)
LRD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $STEP_PARAMS | cut -f5)
TUNE_EVERY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $STEP_PARAMS | cut -f6)

# make outdir 

OUTDIR=training_ALT_Ep3

if [ ! -d "$WD/${OUTDIR}" ]; then
   mkdir -p $WD/${OUTDIR}/
fi

MODEL_SUBDIR=/training_ALT/checkpoints/ckpt-750

## Initial model to start training from

apptainer run \
--nv \
-B $WD:/home \
$DV_PATH \
/opt/deepvariant/bin/train \
--config=/home/dv_config.py:base \
--config.train_dataset_pbtxt="/home/examples_shuffled_ALT/train/shuf_3/examples_shuf3_config.pbtxt" \
--config.tune_dataset_pbtxt="/home/examples_shuffled_ALT/tune/tune_examples.config.pbtxt" \
--config.num_epochs=2 \
--config.learning_rate=$LR \
--config.learning_rate_decay_rate=$LRD \
--config.num_validation_examples=0 \
--config.tune_every_steps=250 \
--experiment_dir=/home/${OUTDIR} \
--strategy=mirrored \
--config.batch_size=$BS \
--config.init_checkpoint="/home/${MODEL_SUBDIR}"


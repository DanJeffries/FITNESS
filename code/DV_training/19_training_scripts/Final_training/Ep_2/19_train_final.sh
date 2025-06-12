#!/bin/bash

#SBATCH --partition=gpu
#SBATCH --time=12:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=2
##SBATCH --gres=gpu:h100:1
#SBATCH --gres=gpu:rtx4090:1
##SBATCH --gres=gpu:rtx3090:1
#SBATCH --mem-per-cpu=40G
#SBATCH --export=NONE
#SBATCH --array=5  # start from 2 because of the header in the params file. 
#SBATCH --job-name=TRAIN_FINAL_EP2
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

####script to run make_examples
WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/

DV_PATH=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV/deepvariant_1.6.1-gpu.modified.sif
OPENBLAS_NUM_THREADS=1 #Set number of threads that OPENBLAS uses to avoid thread overflow error in numpy

# get step instructions

STEP_PARAMS="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/19_training_scripts/Final_training/Ep_2/final_parameters.txt"

STEP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $STEP_PARAMS | cut -f1)
RUN=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $STEP_PARAMS | cut -f2)
BS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $STEP_PARAMS | cut -f3)
LR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $STEP_PARAMS | cut -f4)
LRD=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $STEP_PARAMS | cut -f5)
TUNE_EVERY=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $STEP_PARAMS | cut -f6)

## get the best checkpoint from the previous epoch

CHECKPOINTS="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/19_training_scripts/Final_training/Ep_2/best_checkpoints_from_Ep1.txt"

CHECKPOINT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $CHECKPOINTS | cut -f2)

export APPTAINER_TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/apptainer_tmp_${STEP}_${RUN}_$CHECKPOINT" #Set Singularity temporary dir
export TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/parralel_tmp_${STEP}_${RUN}_$CHECKPOINT" #Set global temporary dir for parallel

# make outdir 

OUTDIR=FINAL_TRAINING/RUN_${RUN}/Epoch_2

if [ ! -d "$WD/${OUTDIR}" ]; then
   mkdir -p $WD/${OUTDIR}/
fi

apptainer run \
--nv \
-B $WD:/home \
$DV_PATH \
/opt/deepvariant/bin/train \
--config=/home/dv_config.py:base \
--config.train_dataset_pbtxt="/home/examples_shuffled_NEW/train/training_subset/examples_shuf2_config.pbtxt" \
--config.tune_dataset_pbtxt="/home/examples_shuffled_NEW/tune/tuning_subset/tune_examples.config.pbtxt" \
--config.num_epochs=1 \
--config.learning_rate=$LR \
--config.learning_rate_decay_rate=$LRD \
--config.learning_rate_num_epochs_per_decay=0.15 \
--config.num_validation_examples=0 \
--config.tune_every_steps=$TUNE_EVERY \
--experiment_dir=/home/$OUTDIR \
--strategy=mirrored \
--config.batch_size=${BS} \
--config.init_checkpoint="/home/FINAL_TRAINING/RUN_${RUN}/Epoch_1/checkpoints/${CHECKPOINT}"

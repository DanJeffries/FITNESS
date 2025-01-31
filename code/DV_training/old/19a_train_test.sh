#!/bin/bash

#SBATCH --partition=gpu
#SBATCH --time=24:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:rtx4090:1
#SBATCH --mem-per-cpu=50G
#SBATCH --export=NONE
##SBATCH --array=1-4
#SBATCH --job-name=TRAIN_tests_learning_rate_decay
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

####script to run make_examples
WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/
PARAMS=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/learning_rates.txt

#LR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $PARAMS | cut -f1)
#BS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $PARAMS | cut -f2)

echo "LR = ${LR}"

LR=0.01
BS=512

export APPTAINER_TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/apptained_tmp/test_${LR}_${BS}" # Set Singularity temporary dir
export TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/parralel_tmp/test_${LR}_${BS}" # Set global temporary dir for parallel

if [ ! -d "$APPTAINER_TMPDIR" ]; then
   mkdir -p $APPTAINER_TMPDIR
fi

if [ ! -d "$TMPDIR" ]; then
   mkdir -p $TMPDIR
fi


DV_PATH=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV/deepvariant_1.6.1-gpu.modified.sif
OPENBLAS_NUM_THREADS=1 # Set number of threads that OPENBLAS uses to avoid thread overflow error in numpy

# make output dir

OUTDIR=training_outs/param_tests/test_${LR}_${BS}/

if [ ! -d "$WD/${OUTDIR}" ]; then
   mkdir -p $WD/${OUTDIR}/
fi

apptainer run \
--nv \
-B $WD:/home \
$DV_PATH \
/opt/deepvariant/bin/train \
--config=/home/dv_config.py:base \
--config.train_dataset_pbtxt="/home/examples_shuffled/train/shuf_test/examples_shuf3_testset_config.pbtxt" \
--config.tune_dataset_pbtxt="/home/examples_shuffled/tune_test/tune_test_examples_config.pbtxt" \
--config.num_epochs=3 \
--config.learning_rate=${LR} \
--config.num_validation_examples=0 \
--config.tune_every_steps=200 \
--config.warmup_steps=200 \
--config.learning_rate_decay_rate=0.9 \
--config.learning_rate_num_epochs_per_decay=1 \
--experiment_dir=/home/${OUTDIR} \
--strategy=mirrored \
--config.batch_size=${BS} \
--config.init_checkpoint="/home/model_wgs_v1.6.1/deepvariant.wgs.ckpt"



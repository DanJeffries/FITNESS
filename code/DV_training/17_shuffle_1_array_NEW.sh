#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=04:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G
#SBATCH --export=NONE
#SBATCH --array=1-10
#SBATCH --job-name=SHUFFLE_train_1_NEW
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load Anaconda3

## Ok this shuffle step is a pain in the arse! Takes so much memory. So I have to split it up and shuffle in two stages. I will split the 200 examples files into 10 batches where each one has a random combo of the files. 

## I will move each batch into a sub-folder and shuffle there. 

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training
SHUFFLE_SCRIPT_DIR=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/UTILS

EXAMPLES_DIR=${WD}/examples_NEW/train
OUTPUT_DIR=${WD}/examples_shuffled_NEW/train/shuf_1

if [ ! -d "$OUTPUT_DIR" ]; then
   mkdir -p $OUTPUT_DIR
fi

## Make the subset files and directories

EXAMPLE_SUBSET_DIR=${WD}/examples_ALT/train/subset_${SLURM_ARRAY_TASK_ID}

## activate python venv
cd /storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/UTILS
. beam/bin/activate

python3 ${SHUFFLE_SCRIPT_DIR}/shuffle_tfrecords_beam.py \
  --project="FITNESS_TRAIN" \
  --input_pattern_list="${EXAMPLE_SUBSET_DIR}/*tfrecord*00020" \
  --output_pattern_prefix="${OUTPUT_DIR}/examples_shuffle1_${SLURM_ARRAY_TASK_ID}" \
  --output_dataset_name="Shuffle_global" \
  --output_dataset_config_pbtxt="${OUTPUT_DIR}/examples_shuffle1_${SLURM_ARRAY_TASK_ID}_config.pbtxt" \
  --job_name=shuffle-tfrecords \
  --runner=DirectRunner \
  --direct_num_workers=20


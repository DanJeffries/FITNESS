#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=3-00:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=990G
#SBATCH --export=NONE
#SBATCH --job-name=SHUFFLE_train_test_large
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load Anaconda3

## Here I will make a test dataset for experimenting with different hyperparameters. 
## Two subsets take 24 hours, so I need to be using a 24th of that really, so its 1 hour per run ideally. 

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training
SHUFFLE_SCRIPT_DIR=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/UTILS

EXAMPLES_DIR=${WD}/examples_shuffled/train/shuf3 ## starting from outputs from shuffle 3.
OUTPUT_DIR=${WD}/examples_shuffled/train/shuf_test_large

if [ ! -d "$OUTPUT_DIR" ]; then
   mkdir -p $OUTPUT_DIR
fi

## Make the subset files and directories

SHUFFLED_FILE_LIST=${WD}/examples_shuffled/train/shuf1/files_shuf.txt

## activate python venv
cd /storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/UTILS
. beam/bin/activate

## making a subset using only shards 00000 - 00010 of 20, i.e. 10 shards, to reduce number of examples. should equate to about 750K examples and take maybe 15 hrs per job. 

python3 ${SHUFFLE_SCRIPT_DIR}/shuffle_tfrecords_beam.py \
  --project="FITNESS_TRAIN" \
  --input_pattern_list="${EXAMPLES_DIR}/examples_shuf3_sub_1_10.shuffled-0000*-of-00020.tfrecord.gz" \
  --output_pattern_prefix="${OUTPUT_DIR}/examples_shuf3_largetestset.shuffled" \
  --output_dataset_name="Shuffle_global" \
  --output_dataset_config_pbtxt="${OUTPUT_DIR}/examples_shuf3_largetestset_config.pbtxt" \
  --job_name=shuffle-tfrecords \
  --runner=DirectRunner \
  --direct_num_workers=20


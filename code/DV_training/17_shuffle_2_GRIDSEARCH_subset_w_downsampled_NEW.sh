#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=12:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=900G
#SBATCH --export=NONE
#SBATCH --job-name=SHUFFLE_full_train_2_subset
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load Anaconda3

## Ok so up to here there have been 2 shuffle steps, so the loci should be pretty well shuffled now. 
## I actually don't need / can't use the whole example dataset that we have created. Training just takes too long. So I figure that what I will do is create a subset via a third shuffle step. This has the advantage of also outputting the appropriate pbtxt file. 

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training
SHUFFLE_SCRIPT_DIR=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/UTILS

EXAMPLES_DIR=${WD}/examples_shuffled_NEW/train/GRIDSEARCH_training_subset_with_downsampled ## starting from outputs from shuffle 1 this time, no shuffle 2 for this test.
OUTPUT_DIR=${WD}/examples_shuffled_NEW/train/GRIDSEARCH_training_subset_with_downsampled_shuffled

if [ ! -d "$OUTPUT_DIR" ]; then
   mkdir -p $OUTPUT_DIR
fi

## activate python venv
cd /storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/UTILS
. beam/bin/activate

python3 ${SHUFFLE_SCRIPT_DIR}/shuffle_tfrecords_beam.py \
  --project="FITNESS_TRAIN" \
  --input_pattern_list="${EXAMPLES_DIR}/examples*tfrecord.gz" \
  --output_pattern_prefix="${OUTPUT_DIR}/examples_gridsearch_with_downsampled_shuffled" \
  --output_dataset_name="Shuffle_global" \
  --output_dataset_config_pbtxt="${OUTPUT_DIR}/examples_gridsearch_with_downsampled_shuffled.pbtxt" \
  --job_name=shuffle-tfrecords \
  --runner=DirectRunner \
  --direct_num_workers=20


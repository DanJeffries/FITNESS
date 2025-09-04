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

EXAMPLES_DIR=${WD}/examples_shuffled_NEW/train/shuf_1 ## starting from outputs from shuffle 1 this time, no shuffle 2 for this test.
OUTPUT_DIR=${WD}/examples_shuffled_NEW/train/GRIDSEARCH_training_subset

if [ ! -d "$OUTPUT_DIR" ]; then
   mkdir -p $OUTPUT_DIR
fi

## Make the subset files and directories

## I first make a randomised list of all of the example files using. Note I can't do this in the array jobs as it would get remade each time, and the random order would be different each time. That would lead to some example files being used twice, and some not at all.

#> ls $EXAMPLES_DIR/*20 | shuf > $EXAMPLES_DIR/shuffled_examples.txt  ## only want file names, no need for full paths

SHUFFLED_EXAMPLES_FILE=$EXAMPLES_DIR/shuffled_examples.txt

## Now we can take every Nth line of this file to create batches. I will split into 10 batches here. MAKE SURE THIS NUMBER IS THE SAME AS THE ARRAY SIZE

## move each file in the subset to its own directory

EXAMPLE_SUBSET_DIR=$EXAMPLES_DIR/gridsearch_subset

if [ ! -d "$EXAMPLE_SUBSET_DIR" ]; then
   mkdir -p $EXAMPLE_SUBSET_DIR
fi

## Note - need the * wild card at the end so the json file gets moved too
for i in $(cat $SHUFFLED_EXAMPLES_FILE | awk "NR % 10 == 2"); do mv ${EXAMPLES_DIR}/${i}* $EXAMPLE_SUBSET_DIR; done

## activate python venv
cd /storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/UTILS
. beam/bin/activate

## making a subset using shuffle subsets 1 and 10 (just for ease of use with the wildcard below)

python3 ${SHUFFLE_SCRIPT_DIR}/shuffle_tfrecords_beam.py \
  --project="FITNESS_TRAIN" \
  --input_pattern_list="${EXAMPLE_SUBSET_DIR}/examples*tfrecord.gz" \
  --output_pattern_prefix="${OUTPUT_DIR}/examples_shuf2_gridsearch" \
  --output_dataset_name="Shuffle_global" \
  --output_dataset_config_pbtxt="${OUTPUT_DIR}/examples_shuf2_gridsearch_config.pbtxt" \
  --job_name=shuffle-tfrecords \
  --runner=DirectRunner \
  --direct_num_workers=20


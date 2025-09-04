#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=06:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=999G
#SBATCH --export=NONE
#SBATCH --array=1-10
#SBATCH --job-name=SHUFFLE_train_1_NEW
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load Anaconda3

## Ok this shuffle step is a pain in the arse! Takes so much memory. So I have to split it up and shuffle in two stages. I will split the 100 examples files into 10 batches where each one has a random combo of the files. 

## I will move each batch into a sub-folder and shuffle there. 

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training
SHUFFLE_SCRIPT_DIR=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/UTILS

EXAMPLES_DIR=${WD}/examples_NEW/train
OUTPUT_DIR=${WD}/examples_shuffled_NEW/train/shuf_1

if [ ! -d "$OUTPUT_DIR" ]; then
   mkdir -p $OUTPUT_DIR
fi


## I first make a randomised list of all of the example files using. Note I can't do this in the array jobs as it would get remade each time, and the random order would be different each time. That would lead to some example files being used twice, and some not at all. 

#> ls $EXAMPLES_DIR/*20 | shuf > $EXAMPLES_DIR/shuffled_examples.txt  ## only want file names, no need for full paths

SHUFFLED_EXAMPLES_FILE=$EXAMPLES_DIR/shuffled_examples.txt

## Now we can take every Nth line of this file to create batches. I will split into 10 batches here. MAKE SURE THIS NUMBER IS THE SAME AS THE ARRAY SIZE 

SUBSET_NUMBER=$(($SLURM_ARRAY_TASK_ID-1)) ## the subset number needs to start from 0
EXAMPLE_SUBSET_DIR=$EXAMPLES_DIR/subset_${SLURM_ARRAY_TASK_ID}

if [ ! -d "$EXAMPLE_SUBSET_DIR" ]; then
   mkdir -p $EXAMPLE_SUBSET_DIR
fi

## move each file in the subset to its own directory. 
## Note - need the * wild card at the end so the json file gets moved too
for i in $(cat $SHUFFLED_EXAMPLES_FILE | awk "NR % 10 == $SUBSET_NUMBER"); do mv ${EXAMPLES_DIR}/${i}* $EXAMPLE_SUBSET_DIR; done

## Make the subset files and directories

EXAMPLE_SUBSET_DIR=${WD}/examples_NEW/train/subset_${SLURM_ARRAY_TASK_ID}

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


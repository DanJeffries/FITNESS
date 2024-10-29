#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=3-00:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=990G
#SBATCH --export=NONE
#SBATCH --array=1-10
#SBATCH --job-name=SHUFFLE_train_1_pt1
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load Anaconda3

## Ok this shuffle step is a pain in the arse! Takes so much memory. So I have to split it up and shuffle in two stages. I will split the 200 examples files into 10 batches where each one has a random combo of the files. 

## I will move each batch into a sub-folder and shuffle there. 

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training
SHUFFLE_SCRIPT_DIR=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/UTILS

EXAMPLES_DIR=${WD}/examples/train
OUTPUT_DIR=${WD}/examples_shuffled/train

if [ ! -d "$OUTPUT_DIR" ]; then
   mkdir -p $OUTPUT_DIR
fi

## Make the subset files and directories

EXAMPLE_SUBSET_DIR=${WD}/examples/train/subset_${SLURM_ARRAY_TASK_ID}

if [ ! -d "$EXAMPLE_SUBSET_DIR" ]; then
   mkdir -p $EXAMPLE_SUBSET_DIR
fi

SHUFFLED_FILE_LIST=${WD}/examples/train/files_shuf.txt

FILE_SUBSET=$(sed -n "${SLURM_ARRAY_TASK_ID}~10p" $SHUFFLED_FILE_LIST)

#for i in $FILE_SUBSET
#do
#	mv ${WD}/examples/train/${i}* $EXAMPLE_SUBSET_DIR/ ## move the record and the json file to new dir
#done

## activate python venv
cd /storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/UTILS
. beam/bin/activate


python3 ${SHUFFLE_SCRIPT_DIR}/shuffle_tfrecords_beam.py \
  --project="FITNESS_TRAIN" \
  --input_pattern_list="${EXAMPLE_SUBSET_DIR}/*tfrecord*00020" \
  --output_pattern_prefix="${OUTPUT_DIR}/examples_shuf1_sub_${SLURM_ARRAY_TASK_ID}.shuffled" \
  --output_dataset_name="Shuffle_global" \
  --output_dataset_config_pbtxt="${OUTPUT_DIR}/examples_shuf1_sub_${SLURM_ARRAY_TASK_ID}_config.pbtxt" \
  --job_name=shuffle-tfrecords \
  --runner=DirectRunner \
  --direct_num_workers=20


#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=2-00:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=990G
#SBATCH --export=NONE
#SBATCH --array=18
#SBATCH --job-name=SHUFFLE_train_2_pt1
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load Anaconda3

## Ok this shuffle step is a pain in the arse! Takes so much memory. So I have to split it up and shuffle in two stages. I will split the 100 examples files into 20 batches where each one has a random combo of 5 files. 

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training
SHUFFLE_SCRIPT_DIR=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/UTILS

EXAMPLES_SHUF_1_DIR=${WD}/examples_shuffled/train/shuf_1/
SHUFFLED_FILE_LIST=${EXAMPLES_SHUF_1_DIR}/files_shuf.txt
EXAMPLE_SUBSET_DIR=${EXAMPLES_SHUF_1_DIR}/subset_${SLURM_ARRAY_TASK_ID}

OUTPUT_DIR=${WD}/examples_shuffled/train/shuf_2

if [ ! -d "$OUTPUT_DIR" ]; then
   mkdir -p $OUTPUT_DIR
fi

## activate python venv
cd /storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/UTILS
. beam/bin/activate

python3 ${SHUFFLE_SCRIPT_DIR}/shuffle_tfrecords_beam.py \
  --project="FITNESS_TRAIN" \
  --input_pattern_list="${EXAMPLE_SUBSET_DIR}/*tfrecord.gz" \
  --output_pattern_prefix="${OUTPUT_DIR}/examples_shuffle2_${SLURM_ARRAY_TASK_ID}" \
  --output_dataset_name="Shuffle_global" \
  --output_dataset_config_pbtxt="${OUTPUT_DIR}/examples_shuffle2_${SLURM_ARRAY_TASK_ID}_config.pbtxt" \
  --job_name=shuffle-tfrecords \
  --runner=DirectRunner \
  --direct_num_workers=20


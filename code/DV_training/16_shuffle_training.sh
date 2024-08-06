#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=10:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=250G
#SBATCH --export=NONE
#SBATCH --job-name=SHUFFLE_global
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/
SHUFFLE_SCRIPT_DIR=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/UTILS/

CROSSES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/crossIDs.txt
CROSS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $CROSSES)

EXAMPLES_DIR=${WD}/examples
OUTPUT_DIR=${WD}/examples_shuffled

if [ ! -d "$OUTPUT_DIR" ]; then
   mkdir $OUTPUT_DIR
fi

## activate pynthon venv
cd /storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/UTILS
. beam/bin/activate

python3 ${SHUFFLE_SCRIPT_DIR}/shuffle_tfrecords_beam.py \
  --project="FITNESS_TRAIN" \
  --input_pattern_list="${EXAMPLES_DIR}/train_all/*_training_examples.tfrecord-?????-of-00020" \
  --output_pattern_prefix="${OUTPUT_DIR}/All_samples_training_examples.shuffled" \
  --output_dataset_name="Shuffle_global" \
  --output_dataset_config_pbtxt="${OUTPUT_DIR}/All_samples_training_examples.dataset_config.pbtxt" \
  --job_name=shuffle-tfrecords \
  --runner=DirectRunner \
  --direct_num_workers=20


#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=01:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=40G
#SBATCH --export=NONE
#SBATCH --job-name=SHOW_EXAMPLES
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

####script to run make_examples
WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training
REF=/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_formatted_shortnames.fna

REGIONS_BED=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/training_regions/${CROSS}_train_partitions.bed

# Make temp dirs to be used instead of in /tmp. (need to be in $HOME)

export APPTAINER_TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/apptained_tmp/$SAMPLE" #Set Singularity temporary dir
export TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/parralel_tmp/$SAMPLE" #Set global temporary dir for parallel

if [ ! -d "$APPTAINER_TMPDIR" ]; then
   mkdir -p $APPTAINER_TMPDIR
fi

if [ ! -d "$TMPDIR" ]; then
   mkdir -p $TMPDIR
fi


DV_PATH=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV/deepvariant_1.6.1-gpu.modified.sif
OPENBLAS_NUM_THREADS=1 #Set number of threads that OPENBLAS uses to avoid thread overflow error in numpy

# make output dir

SAMPLE=FG_male_1
REGION=NC_053222.1_chromosome_11:7000-7500

if [ ! -d "$WD/examples/show/${SAMPLE}_${REGION}" ]; then
   mkdir -p $WD/examples/show/${SAMPLE}_${REGION}
fi

apptainer run \
-B $WD:/wd \
$DV_PATH  \
  /opt/deepvariant/bin/show_examples \
  --examples=/wd/examples/train/${SAMPLE}_training_examples_positional.tfrecord@20 \
  --example_info_json=/wd/examples/train/${SAMPLE}_training_examples_positional.tfrecord-00000-of-00020.example_info.json \
  --output=/wd/examples/show/${SAMPLE}_${REGION}/${SAMPLE}_pileup \
  --num_records=20 \
  --curate \
  --regions $REGION \




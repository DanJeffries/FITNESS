#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=02:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G
#SBATCH --export=NONE
#SBATCH --array=1-5
#SBATCH --job-name=MAKE_EX_TRAIN_NEW_DOWNSAMPLE
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

####script to run make_examples
WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/
REF=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_formatted_shortnames.fna
SAMPLES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/offspring.txt
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $SAMPLES)

CROSS=$(echo $SAMPLE | cut -f1 -d'_') ## get cross ID so I can get the regions for this cross

## I created the training regions files manually - note I cant use this bash var in the docker, this is just to show where it is. 

REGIONS_BED=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/training_regions/${CROSS}_train_partitions.bed

# Make temp dirs to be used instead of in /tmp. (need to be in $HOME)

export APPTAINER_TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/apptained_tmp" #Set Singularity temporary dir
export TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/parralel_tmp" #Set global temporary dir for parallel

DV_PATH=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV/deepvariant_1.6.1-gpu.modified.sif
OPENBLAS_NUM_THREADS=1 #Set number of threads that OPENBLAS uses to avoid thread overflow error in numpy

# make output dir

if [ ! -d "$WD/examples/train" ]; then
   mkdir -p $WD/examples/train
fi

SAMPLE_BED_NAME=$(echo $SAMPLE | cut -f-2 -d'_')

apptainer run \
-B $WD:/wd \
$DV_PATH  \
parallel -q --halt 2 --line-buffer \
/opt/deepvariant/bin/make_examples \
--mode training \
--ref $REF \
--reads /wd/bams/${SAMPLE}.fixmate.coordsorted.bam \
--truth_variants /wd/Filtered_variants/${SAMPLE}.ALL_TRUTH_VARS.vcf.gz \
--confident_regions /wd/Confident_regions/${SAMPLE_BED_NAME}_conf_regions_inc_vars.bed \
--examples /wd/examples/train/${SAMPLE}_training_examples_positional.downsampled_05.tfrecord@20 \
--regions /wd/training_regions/${CROSS}_train_partitions.bed \
--labeler_algorithm=positional_labeler \
--downsample_fraction=0.5 \
--channels "insert_size" \
--task {} ::: `seq 0 19` #split the task into 20 jobs




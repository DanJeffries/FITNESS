#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=02:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=110G
#SBATCH --export=NONE
#SBATCH --array=11 #-12
#SBATCH --job-name=MAKE_EX_CHR1_TEST
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

####script to run make_examples
WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/
REF=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_formatted_shortnames.fna
SAMPLES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/parents.txt

SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $SAMPLES)

# Make temp dirs to be used instead of in /tmp.
#export APPTAINER_TMPDIR="/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/apptained_tmp" #Set Singularity temporary dir
#export TMPDIR="/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/parralel_tmp" #Set global temporary dir for parallel

export APPTAINER_TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/apptained_tmp" #Set Singularity temporary dir
export TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/parralel_tmp" #Set global temporary dir for parallel

DV_PATH=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV/deepvariant_1.6.1-gpu.modified.sif
OPENBLAS_NUM_THREADS=1 #Set number of threads that OPENBLAS uses to avoid thread overflow error in numpy

# make output dir

if [ ! -d "$WD/examples/${SAMPLE}" ]; then
   mkdir $WD/examples/${SAMPLE}
fi

apptainer run \
-B $WD:/wd \
$DV_PATH  \
parallel -q --halt 2 --line-buffer \
/opt/deepvariant/bin/make_examples \
--mode training \
--ref $REF \
--reads /wd/bams/${SAMPLE}.fixmate.coordsorted.bam \
--truth_variants /wd/Filtered_variants/${SAMPLE}.vcf.gz \
--confident_regions /wd/Confident_regions/${SAMPLE}.conf.1n_repeats_removed.bed \
--examples /wd/examples/${SAMPLE}/make_examples.tfrecord@20 \
--channels "insert_size" \
--task {} ::: `seq 0 19` #split the task into 20 jobs


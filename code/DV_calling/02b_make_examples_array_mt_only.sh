#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=01:30:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=2G
#SBATCH --export=NONE
#SBATCH --array=1-100
#SBATCH --job-name=MAKE_EXAMPLES_CALL_mt_cohort_1
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

####script to run make_examples
## here i am calling only on the mitochondria 

WD=/storage/research/iee_temp_dj20y461/DV_calling

REF=/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_formatted_shortnames.fna
ITER_FILE=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/sample_metadata/sample_paths_rmdup.txt
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE | rev | cut -f1 -d'/' | rev)
BAM_DIR=bams
BAM_PATH=$BAM_DIR/${SAMPLE}.fixmate.coordsorted.bam 


MASK='mask/exclude_nuclear.bed'

# Make temp dirs to be used instead of in /tmp. (need to be in $SCRATCH for calling as it makes a lot of temp files)

APPTAINER_TMPDIR="$WD/apptained_tmp/${SAMPLE}_calling" #Set Singularity temporary dir

TMPDIR="/parralel_tmp/${SAMPLE}_calling" #Set global temporary dir for parallel - relativ to $WD which will be mounted in the container


if [ ! -d "$APPTAINER_TMPDIR" ]; then
     	mkdir -p $APPTAINER_TMPDIR
fi

if [ ! -d "$WD/$TMPDIR" ]; then
     	mkdir -p $WD/$TMPDIR
fi


DV_PATH=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV/deepvariant_1.6.1-gpu.modified.sif
OPENBLAS_NUM_THREADS=1 #Set number of threads that OPENBLAS uses to avoid thread overflow error in numpy

# make example output dir

if [ ! -d "$WD/examples_mt/$SAMPLE" ]; then
   mkdir -p $WD/examples_mt/$SAMPLE
fi

## make gvcf TF records output dir

if [ ! -d "$WD/GVCF_TFRECORDS_mt/$SAMPLE" ]; then
   mkdir -p $WD/GVCF_TFRECORDS_mt/$SAMPLE
fi

N_SHARDS=20

apptainer run \
-B $WD:/wd \
$DV_PATH  \
parallel -q --halt 2 --line-buffer --tmpdir /wd/$TMPDIR \
/opt/deepvariant/bin/make_examples \
--mode calling \
--ref /wd/$REF \
--reads /wd/$BAM_PATH \
--examples /wd/examples_mt/$SAMPLE/${SAMPLE}_examples_positional.tfrecord@${N_SHARDS}.gz \
--gvcf /wd/GVCF_TFRECORDS_mt/$SAMPLE/${SAMPLE}_gvcf.tfrecord@${N_SHARDS}.gz \
--exclude_regions /wd/$MASK \
--labeler_algorithm=positional_labeler \
--channels "insert_size" \
--task {} ::: `seq 0 19` #split the task into 20 jobs

rm -rf $APPTAINER_TMPDIR*
rm -rf $WD/$TMPDIR*


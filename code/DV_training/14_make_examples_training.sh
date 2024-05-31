#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=02:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=200G
#SBATCH --export=NONE
#SBATCH --array=1-10
#SBATCH --job-name=MAKE_EX_CHR1_TEST
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

####script to run make_examples
WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/
REF=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_formatted_shortnames.fna
SAMPLES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/offspring.txt
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $SAMPLES)

## Get the regions to use for this step

CROSS=$(echo $SAMPLE | cut -f1 -d'_') ## get cross ID so I can get the regions for this cross

## I created the training regions files manually 

REGIONS_BED=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/training_regions/training_regions_$CROSS.bed

for CHROM in $(echo $REGIONS); 
  do
    CHROM_NAME=$(grep $CHROM $FAI | cut -f1)
    END=$(grep $CHROM $FAI | cut -f2)
    echo -e "$CHROM_NAME\t0\t$END"
  
  done > $REGIONS_BED


# Make temp dirs to be used instead of in /tmp. (need to be in $HOME)

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
--truth_variants /wd/Filtered_variants/${SAMPLE}.ALL_TRUTH_VARS.CORRECTED.vcf.gz \
--confident_regions /wd/Confident_regions/${SAMPLE}.conf.bed \
--examples /wd/examples/${SAMPLE}/training_examples.tfrecord@20 \
--regions /wd/training_regions/training_regions_$CROSS.bed \
--channels "insert_size" \
--task {} ::: `seq 0 19` #split the task into 20 jobs


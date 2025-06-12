#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=06:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=2G
#SBATCH --export=NONE
#SBATCH --array=1-18
#SBATCH --job-name=CALL_AND_PROCESS_VARIANTS_filter_calibration
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

####script to call variants from examples

#WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_calling

WD=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_calling/filter_calibration
REF=/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames.fna

## iterating over the same master sample file as for all the other steps
ITER_FILE=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_calling/filter_calibration/family_samples.txt
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE | rev | cut -f1 -d'/' | rev)

EXAMPLE_SHARDS=20
PROCESS_SHARDS=16

EXAMPLES_DIR=examples/$SAMPLE

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

DEEPSTICKLEv1=DeepStickle_v1.0/ckpt-32068

# make call_variants output dir

CALL_VARIANTS_DIR=call_variants/$SAMPLE

if [ ! -d "$WD/$CALL_VARIANTS_DIR" ]; then
   mkdir -p $WD/$CALL_VARIANTS_DIR
fi

#apptainer run \
#-B $WD:/wd \
#$DV_PATH  \
#/opt/deepvariant/bin/call_variants \
#	--checkpoint /wd/${DEEPSTICKLEv1}  \
#	--examples /wd/$EXAMPLES_DIR/${SAMPLE}_examples_positional.tfrecord@${EXAMPLE_SHARDS}.gz  \
#	--outfile /wd/$CALL_VARIANTS_DIR/${SAMPLE}_call_variants.tfrecord.gz \
#	--writer_threads $PROCESS_SHARDS \
#	--num_readers $EXAMPLE_SHARDS \
#
#echo "call variants done"

COHORT=families  ## output to cohort specific subdirs so that I can more easily give the path to GLnexus
POSTPROCESS_DIR=postprocess_variants/$COHORT/$SAMPLE

# make output dir

if [ ! -d "$WD/$POSTPROCESS_DIR/" ]; then
   mkdir -p $WD/$POSTPROCESS_DIR/
fi

apptainer run \
-B $WD:/wd \
$DV_PATH  \
/opt/deepvariant/bin/postprocess_variants \
        --ref "/wd/$REF" \
        --infile "/wd/$CALL_VARIANTS_DIR/${SAMPLE}_call_variants@${PROCESS_SHARDS}.tfrecord.gz" \
        --outfile "/wd/$POSTPROCESS_DIR/${SAMPLE}.vcf.gz" \
	--nonvariant_site_tfrecord_path "/wd/GVCF_TFRECORDS/$SAMPLE/${SAMPLE}_gvcf.tfrecord@${EXAMPLE_SHARDS}.gz" \
	--gvcf_outfile "/wd/$POSTPROCESS_DIR/${SAMPLE}.g.vcf.gz" \
        --cpus "16"

echo "post processing done"

rm -rf $APPTAINER_TMPDIR*
rm -rf $WD/$TMPDIR*


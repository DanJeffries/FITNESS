#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=10:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=5G
##SBATCH --mem=100G
#SBATCH --export=NONE
#SBATCH --array=2 #-100
#SBATCH --job-name=CALL_VARIANTS
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

####script to call variants from examples

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_calling
REF=/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_formatted_shortnames.fna

## dont need bam files now but can use this as the list to get sample names from
BAM_FILES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_calling/bam_batch.txt
BAM=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $BAM_FILES)
SAMPLE=$(echo $BAM | cut -f1 -d'.') 

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

# make output dir

OUTDIR=call_variants

if [ ! -d "$WD/$OUTDIR/$SAMPLE" ]; then
   mkdir -p $WD/$OUTDIR/$SAMPLE
fi

apptainer run \
-B $WD:/wd \
$DV_PATH  \
/opt/deepvariant/bin/call_variants \
	--checkpoint /wd/${DEEPSTICKLEv1}  \
	--examples /wd/examples/$SAMPLE/${SAMPLE}_examples_positional.tfrecord@20.gz  \
	--outfile /wd/$OUTDIR/$SAMPLE/${SAMPLE}_call_variants.tfrecord.gz \
	--writer_threads 16 \
	--num_readers 16 \

#--outfile "/wd/parralel_tmp/FG_CC_19T_011_calling/tmpxvh4ull2/call_variants_output.tfrecord.gz" --examples "/wd/parralel_tmp/FG_CC_19T_011_calling/tmpxvh4ull2/make_examples.tfrecord@20.gz" --checkpoint "/wd/DeepStickle_v1.0/ckpt-32068"


rm -rf $APPTAINER_TMPDIR*
rm -rf $WD/$TMPDIR*


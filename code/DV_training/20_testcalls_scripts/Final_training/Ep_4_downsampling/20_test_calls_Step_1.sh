#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=06:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
##SBATCH --gres=gpu:rtx4090:1
#SBATCH --mem-per-cpu=10G
#SBATCH --export=NONE
#SBATCH --job-name=FINAL_TRAINING_test_calls_Ep4
#SBATCH --array=2-3
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

####script to run make_examples
WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/
REF=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_formatted_shortnames.fna
SAMPLES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/sample_names.txt

SAMPLE="FG_male_1"
CROSS="FG"

DV_PATH=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV/deepvariant_1.6.1-gpu.modified.sif
OPENBLAS_NUM_THREADS=1 #Set number of threads that OPENBLAS uses to avoid thread overflow error in numpy

## Model to test
### this text file gives path to best checkpoint for each of the training runs being evaluated. Paths start from the path mounted as home in apptainer. 
BEST_MODELS=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/20_testcalls_scripts/Final_training/Ep_4_downsampling/best_checkpoints.txt

MODEL_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $BEST_MODELS | cut -f1)
CHECKPOINT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $BEST_MODELS | cut -f2)
MODEL=${MODEL_DIR}/${CHECKPOINT}

## make outdir

MODEL_OUTDIR=${MODEL_DIR}/test_calls

if [ ! -d "$WD/$MODEL_OUTDIR" ]; then
   mkdir -p "$WD/$MODEL_OUTDIR"
fi

## set temp dirs for the run
export APPTAINER_TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/apptained_tmp/model_tests_${CHECKPOINT}_${SAMPLE}" #Set Singularity temporary dir
export TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/parralel_tmp/model_tests_${CHECKPOINT}_${SAMPLE}" #Set global temporary dir for parallel

if [ ! -d "$APPTAINER_TMPDIR" ]; then
   mkdir -p "$APPTAINER_TMPDIR"
fi

if [ ! -d "$TMPDIR" ]; then
   mkdir -p "$TMPDIR"
fi

apptainer run \
-B $WD:/home \
$DV_PATH \
/opt/deepvariant/bin/run_deepvariant \
--model_type WGS \
--customized_model "/home/$MODEL" \
--ref "/home/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_formatted_shortnames.fna" \
--reads "/home/bams/${SAMPLE}.fixmate.coordsorted.bam" \
--regions "/home/training_regions/${CROSS}_test_partitions.bed" \
--output_vcf "/home/${MODEL_OUTDIR}/${SAMPLE}_${CHECKPOINT}_test_calls.vcf.gz" \
--num_shards=20


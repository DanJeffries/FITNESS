#!/bin/bash

#SBATCH --partition=gpu
#SBATCH --time=24:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:rtx4090:1
#SBATCH --mem-per-cpu=50G
#SBATCH --export=NONE
#SBATCH --job-name=TRAIN_test1
#SBATCH --array=2-5
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

####script to run make_examples
WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/
REF=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_formatted_shortnames.fna
SAMPLES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/offspring.txt
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $SAMPLES)
CROSS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $SAMPLES | cut -f1 -d'_')

export APPTAINER_TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/apptained_tmp" #Set Singularity temporary dir
export TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/parralel_tmp" #Set global temporary dir for parallel

DV_PATH=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV/deepvariant_1.6.1-gpu.modified.sif
OPENBLAS_NUM_THREADS=1 #Set number of threads that OPENBLAS uses to avoid thread overflow error in numpy

# make output dir

if [ ! -d "$WD/examples/validate/${SAMPLE}" ]; then
   mkdir -p $WD/examples/validate/${SAMPLE}
fi

MODEL_SUBDIR=/training_outs/epoch1/sub12/checkpoints/ckpt-10000

apptainer run \
--nv \
-B $WD:/home \
$DV_PATH \
/opt/deepvariant/bin/run_deepvariant \
--model_type WGS \
--customized_model "/home/${MODEL_SUBDIR}" \
--ref "/home/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_formatted_shortnames.fna" \
--reads "/home/bams/${SAMPLE}.fixmate.coordsorted.bam" \
--regions "/home/training_regions/${CROSS}_test_partitions.bed" \
--output_vcf "/home/test/${SAMPLE}_test_set.vcf.gz" \
--num_shards=20


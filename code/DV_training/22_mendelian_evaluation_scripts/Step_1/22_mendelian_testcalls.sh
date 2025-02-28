#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=06:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
##SBATCH --gres=gpu:rtx4090:1
#SBATCH --mem-per-cpu=10G
#SBATCH --export=NONE
#SBATCH --job-name=MEND_TEST_calls_STEP_1
#SBATCH --array=2-11,13-18 
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

DV_PATH=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV/deepvariant_1.6.1-gpu.modified.sif
OPENBLAS_NUM_THREADS=1 #Set number of threads that OPENBLAS uses to avoid thread overflow error in numpy

####script to run make_examples
WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/
REF=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_formatted_shortnames.fna


## Need to call each sample from the pedigree with each model we are testing. So here I made a model-sample map using the below one-liner. 
# for i in $(sed 's/\t/XXX/g' best_checkpoints.txt) ; do for j in $(cat sample_names.txt) ;do echo $i $j | sed 's/XXX/ /g' ; done; done > model_sample_map.txt

MODEL_SAMPLE_MAP=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/22_mendelian_evaluation_scripts/Step_1/model_sample_map.txt

MODEL_DIR=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $MODEL_SAMPLE_MAP | cut -f1 -d' ')
CHECKPOINT=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $MODEL_SAMPLE_MAP | cut -f2 -d' ')
MODEL=${MODEL_DIR}/${CHECKPOINT}
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $MODEL_SAMPLE_MAP | cut -f3 -d' ')
CROSS=$(echo $SAMPLE | cut -f1 -d'_')

STEP=$(echo $MODEL_DIR | cut -f2 -d"/")
RUN=$(echo $MODEL_DIR | cut -f3 -d"/")

echo "INPUTS ------- "
echo "MODEL_DIR $MODEL_DIR"
echo "CHECKPOINT $CHECKPOINT"
echo "MODEL $MODEL"
echo "SAMPLE $SAMPLE"
echo "CROSS $CROSS"
echo "STEP $STEP"
echo "RUN $RUN"


## make outdir

MODEL_OUTDIR=${MODEL_DIR}/test_calls

if [ ! -d "$WD/$MODEL_OUTDIR" ]; then
   mkdir -p "$WD/$MODEL_OUTDIR"
fi

## set temp dirs for the run
export APPTAINER_TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/apptained_tmp/model_tests_${STEP}_${RUN}_${CHECKPOINT}_${SAMPLE}" #Set Singularity temporary dir
export TMPDIR="/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/parralel_tmp/model_tests_${STEP}_${RUN}_${CHECKPOINT}_${SAMPLE}" #Set global temporary dir for parallel

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
--output_vcf "/home/${MODEL_OUTDIR}/${SAMPLE}_${STEP}_${RUN}_${CHECKPOINT}_test_calls.vcf.gz" \
--num_shards=20


#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --export=NONE
#SBATCH --array=2-37
#SBATCH --job-name=FIT_MultiQC
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load MultiQC/1.11-foss-2021a

DATA_DIR=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/raw_BACKUP/
RUN_NAMES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_calling/scripts/runs.txt

RUN_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $RUN_NAMES)
RUN_PATH=${DATA_DIR}/${RUN_NAME}

MULTIQC_RUN_OUTDIR=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/MULTIQC/$RUN_NAME

if [ ! -d "$MULTIQC_RUN_OUTDIR" ]
then
    mkdir -p $MULTIQC_RUN_OUTDIR
fi  

multiqc --ignore *_I[1-2]_* $RUN_PATH -o $MULTIQC_RUN_OUTDIR


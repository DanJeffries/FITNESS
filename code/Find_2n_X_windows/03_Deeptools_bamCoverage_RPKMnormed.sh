#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --export=NONE
#SBATCH --array=1-25
#SBATCH --job-name=Ga_DEEP_500
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module add UHTS/Analysis/deepTools/2.5.4;

ID=$SLURM_ARRAY_TASK_ID

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/Find_2n_X_windows
ITERFILE=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/Find_2n_X_windows/samples_new.txt

SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITERFILE)

BAMDIR=$WD/bams
BAM=${BAMDIR}/${SAMPLE_NAME}.fixmate.coordsorted.bam

OUTDIR=$WD/depths

if [ ! -d "$OUTDIR" ]; then
  mkdir $OUTDIR
fi

bamCoverage --bam $BAM \
            --binSize 500 \
            --numberOfProcessors 16 \
            --verbose \
	    --normalizeUsingRPKM \
            --outFileName ${OUTDIR}/${SAMPLE_NAME}.500bp.RPKM.depth \
            --outFileFormat bedgraph





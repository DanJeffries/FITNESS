#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH --export=NONE
#SBATCH --array=1-49
#SBATCH --job-name=Ga_DEEP_500
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module add UHTS/Analysis/deepTools/2.5.4;

ID=$SLURM_ARRAY_TASK_ID
ITER_FILE=/storage/scratch/iee/dj20y461/Stickleback/Y_comp/Sex_chrom_stats/Read_depths/Ga/scripts/samples.txt
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE)

MARKED_DUPS_DIR=/storage/scratch/iee/dj20y461/Stickleback/Y_comp/Sex_chrom_stats/Read_depths/Ga/bam
MARKED_DUPS_BAM=${MARKED_DUPS_DIR}/${SAMPLE_NAME}.markdup.q20.bam

OUTDIR=/storage/scratch/iee/dj20y461/Stickleback/Y_comp/Sex_chrom_stats/Read_depths/Ga/depths/

if [ ! -d "$OUTDIR" ]; then
  mkdir $OUTDIR
fi


bamCoverage --bam $MARKED_DUPS_BAM \
            --binSize 500 \
            --numberOfProcessors 16 \
            --verbose \
	    --normalizeUsingRPKM \
            --outFileName ${OUTDIR}/${SAMPLE_NAME}.500bp.RPKM.depth \
            --outFileFormat bedgraph





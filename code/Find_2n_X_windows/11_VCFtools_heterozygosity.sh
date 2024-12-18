#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --export=NONE
#SBATCH --job-name=VCFtools_Hardy
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.err

module load vital-it
module add UHTS/Analysis/vcftools/0.1.15;

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/Find_2n_X_windows/
INVCF=$WD/Filtered_VCF/Filtered_BiAL_QUAL_GQ30_MinDP10_AF085_HARD.vcf

MALES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/Find_2n_X_windows/males.txt
FEMALES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/Find_2n_X_windows/females.txt

OUTDIR=$WD/VCF_STATS
OUT_PREFIX=$OUTDIR/Filtered_BiAL_QUAL_GQ30_MinDP10_AF085_HARD

if [ ! -d "$OUTDIR" ]; then
   mkdir $OUTDIR
fi


vcftools \
    --gzvcf $INVCF \
    --keep $MALES \
    --max-missing-count 10 \
    --hardy \
    --out ${OUT_PREFIX}_MALES

vcftools \
    --gzvcf $INVCF \
    --keep $FEMALES \
    --max-missing-count 10 \
    --hardy \
    --out ${OUT_PREFIX}_FEMALES


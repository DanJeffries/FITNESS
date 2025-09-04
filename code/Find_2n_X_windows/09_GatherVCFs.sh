#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --export=NONE
#SBATCH --job-name=Gather
#SBATCH --output=%x_%A.out
#SBATCH --error=%x_%A.err

module load vital-it
module load UHTS/Analysis/GenomeAnalysisTK/4.1.3.0
export GATK_LOCAL_JAR=/software/UHTS/Analysis/GenomeAnalysisTK/4.1.3.0/bin/GenomeAnalysisTK.jar

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/Find_2n_X_windows/
VCFs_LIST=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/Find_2n_X_windows/vcf_paths.list
OUTDIR=$WD/Unfiltered_VCF

if [ ! -d "$OUTDIR" ]; then
   mkdir $OUTDIR
fi

gatk --java-options "-Xmx24g -Xms4g" GatherVcfs \
   -I $VCFs_LIST \
   -O $OUTDIR/Joint_all_unfiltered.vcf



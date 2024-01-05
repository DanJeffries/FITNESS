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

VCFs_PATH=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/
OUTDIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/Unfiltered_VCF/

if [ ! -d "$OUTDIR" ]; then
   mkdir $OUTDIR
fi


gatk --java-options "-Xmx24g -Xms4g" GatherVcfs \
   -I $VCFs_PATH/Joint_called_VCFs_chrI/chrI.vcf \
   -I $VCFs_PATH/Joint_called_VCFs_chrII/chrII.vcf \
   -I $VCFs_PATH/Joint_called_VCFs_chrIII/chrIII.vcf \
   -I $VCFs_PATH/Joint_called_VCFs_chrIV/chrIV.vcf \
   -I $VCFs_PATH/Joint_called_VCFs_chrIX/chrIX.vcf \
   -I $VCFs_PATH/Joint_called_VCFs_chrV/chrV.vcf \
   -I $VCFs_PATH/Joint_called_VCFs_chrVI/chrVI.vcf \
   -I $VCFs_PATH/Joint_called_VCFs_chrVIII/chrVIII.vcf \
   -I $VCFs_PATH/Joint_called_VCFs_chrVII/chrVII.vcf \
   -I $VCFs_PATH/Joint_called_VCFs_chrXI/chrXI.vcf \
   -I $VCFs_PATH/Joint_called_VCFs_chrXII/chrXII.vcf \
   -I $VCFs_PATH/Joint_called_VCFs_chrXIII/chrXIII.vcf \
   -I $VCFs_PATH/Joint_called_VCFs_chrXIV/chrXIV.vcf \
   -I $VCFs_PATH/Joint_called_VCFs_chrXIX/chrXIX.vcf \
   -I $VCFs_PATH/Joint_called_VCFs_chrX/chrX.vcf \
   -I $VCFs_PATH/Joint_called_VCFs_chrXVI/chrXVI.vcf \
   -I $VCFs_PATH/Joint_called_VCFs_chrXVIII/chrXVIII.vcf \
   -I $VCFs_PATH/Joint_called_VCFs_chrXVII/chrXVII.vcf \
   -I $VCFs_PATH/Joint_called_VCFs_chrXV/chrXV.vcf \
   -I $VCFs_PATH/Joint_called_VCFs_chrXX/chrXX.vcf \
   -I $VCFs_PATH/Joint_called_VCFs_chrXXI/chrXXI.vcf \
   -O $OUTDIR/Joint_called.vcf



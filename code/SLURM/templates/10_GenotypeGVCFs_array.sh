#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --job-name=GT_GVCFs
#SBATCH --array=1-21
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Analysis/GenomeAnalysisTK/4.1.3.0
export GATK_LOCAL_JAR=/software/UHTS/Analysis/GenomeAnalysisTK/4.1.3.0/bin/GenomeAnalysisTK.jar

INTERVALS=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/scripts/All_intervals.list
INTERVAL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $INTERVALS)

ASSEMBLY=/storage/scratch/iee/dj20y461/Gacu_assembly/No_Y/stickleback_v5_assembly_NoY.fa
GENDB=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/GenDB_${INTERVAL}
OUTDIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/Joint_called_VCFs_${INTERVAL}/

if [ ! -d $OUTDIR ]; then
   mkdir $OUTDIR
fi


gatk --java-options "-Xmx24g -Xms4g" GenotypeGVCFs \
   -R ${ASSEMBLY} \
   -V gendb://${GENDB} \
   -O ${OUTDIR}/${INTERVAL}.vcf 




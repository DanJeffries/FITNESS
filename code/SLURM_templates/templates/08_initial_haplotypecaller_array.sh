#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=05:00:00
#eBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --array=2-1932
#SBATCH --job-name=HapCall
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Analysis/GenomeAnalysisTK/4.1.3.0

ID=$SLURM_ARRAY_TASK_ID
ITER_FILE=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/scripts/HapCaller_iterfile.txt

SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE | cut -f1 -d' ')
INTERVAL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE | cut -f2 -d' ')

GENOME=/storage/scratch/iee/dj20y461/Gacu_assembly/No_Y/stickleback_v5_assembly_NoY.fa

DATA_DIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/bamMDRG/
OUTDIR=/storage/scratch/iee/dj20y461/FITNESS/MCGILL_TEST_ALL/GVCFs/

if [ ! -d "$OUTDIR" ]; then
   mkdir $OUTDIR
fi

SAMPLE=${DATA_DIR}/${SAMPLE_NAME}.MD.RG.bam

OUTFILE=${OUTDIR}/${SAMPLE_NAME}_${INTERVAL}.g.vcf.gz

java -jar $GATK_PATH/bin/GenomeAnalysisTK.jar HaplotypeCaller \
	-R $GENOME \
	-I $SAMPLE \
	-L $INTERVAL \
	-ERC GVCF \
	--output $OUTFILE



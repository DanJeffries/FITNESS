#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --array=1
#SBATCH --job-name=GenDB_INT_array
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Analysis/GenomeAnalysisTK/4.1.3.0
export GATK_LOCAL_JAR=/software/UHTS/Analysis/GenomeAnalysisTK/4.1.3.0/bin/GenomeAnalysisTK.jar

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/Find_2n_X_windows/

SAMPLE_MAP_TEMPLATE=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/Find_2n_X_windows/gvcf_sample_maps/sample_map_template.txt
INTERVALS=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/Find_2n_X_windows/Intervals.txt
INTERVAL=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $INTERVALS)

SAMPLE_MAP=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/Find_2n_X_windows/gvcf_sample_maps/GVCF_sample_map_${INTERVAL}.txt

sed "s/XXX/$INTERVAL/g" $SAMPLE_MAP_TEMPLATE > $SAMPLE_MAP

GenDB_DIR=$WD/GenDBs/GenDB_${INTERVAL}/

TMP_DIR=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/Find_2n_X_windows/

if [ -d $GenDB_DIR ]; then
   rm ${GenDB_DIR}* -rf
   rmdir $GenDB_DIR
fi

READER_THREADS=16

gatk --java-options "-Xmx128g -Xms128g" GenomicsDBImport \
					--sample-name-map $SAMPLE_MAP \
					--tmp-dir $TMP_DIR \
					--intervals $INTERVAL\
					--genomicsdb-workspace-path $GenDB_DIR \
					--reader-threads $READER_THREADS \
					--consolidate \
					--verbosity DEBUG



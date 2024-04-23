#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=60G
#SBATCH --export=NONE
#SBATCH --array=1-25
#SBATCH --job-name=BWA
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Aligner/bwa/0.7.17
module load UHTS/Analysis/samblaster/0.1.24 
module load UHTS/Analysis/samtools/1.10
module load UHTS/Analysis/sambamba/0.7.1

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/Find_2n_X_windows
ITERFILE=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/Find_2n_X_windows/samples_new.txt

DATA_DIR=$WD/trimmed
SAMPLE_NAME=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITERFILE)

SAMPLE_R1=$WD/trimmed/${SAMPLE_NAME}.R1.trimmed.fastq.gz
SAMPLE_R2=$WD/trimmed/${SAMPLE_NAME}.R2.trimmed.fastq.gz
GENOME_IDX=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/ref/BWA/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY

## Get the read group info to add to bams from the first read of every file

HEADER=$(zcat $SAMPLE_R1 | head -n1) ## get first read ID line 

INSTRUMENT=$(echo $HEADER | cut -f1 -d':'| sed 's/@//g')
RUN=$(echo $HEADER | cut -f2 -d':')
FLOWCELL=$(echo $HEADER | cut -f3 -d':')
LANE=$(echo $HEADER | cut -f4 -d':')
BARCODES=$(echo $HEADER | cut -f2 -d' '| cut -f4 -d':' | sed 's/+/-/g')

echo "@RG\tID:${SAMPLE_NAME}\tPL:ILLUMINA\tPU:${FLOWCELL}.${LANE}\tSM:${SAMPLE_NAME}\tCN:MCGILL\tBC:${BARCODES}"

echo $READGROUP_STRING

OUTDIR=$WD/bams/

if [ ! -d "$OUTDIR" ]; then
   mkdir $OUTDIR
fi

BAM_FIX_COORDSORT=${OUTDIR}/${SAMPLE_NAME}.fixmate.coordsorted.bam

bwa mem -t 20 \
        -R $(echo "@RG\tID:${SAMPLE_NAME}\tPL:ILLUMINA\tPU:${FLOWCELL}.${LANE}\tSM:${SAMPLE_NAME}\tCN:MCGILL\tBC:${BARCODES}") \
        $GENOME_IDX \
        $SAMPLE_R1 \
        $SAMPLE_R2 | \
samblaster | \
samtools fixmate \
	-m \
        -@ 20 \
        - \
	- | \
samtools sort \
	-@ 20 \
	-O BAM \
	-o $BAM_FIX_COORDSORT \
	-

samtools index -c -@ 20 $BAM_FIX_COORDSORT

DEPTHS_OUT=${OUTDIR}/${SAMPLE_NAME}.samcov

samtools coverage $BAM_FIX_COORDSORT > $DEPTHS_OUT

## Calculate some averages and add to the end of each file
egrep -v '#|chrUn|chrM' $DEPTHS_OUT | awk '{ sum += $6; n++ } END { if (n > 0) print "\n## Avg. perc. bases. covered = " sum / n; }' >> $DEPTHS_OUT
egrep -v '#|chrUn|chrM' $DEPTHS_OUT | awk '{ sum += $7; n++ } END { if (n > 0) print "## Avg. coverage = " sum / n; }' >> $DEPTHS_OUT


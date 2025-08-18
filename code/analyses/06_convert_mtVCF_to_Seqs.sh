#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --export=NONE
#SBATCH --array=1-100
#SBATCH --job-name=mt_consensus
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load SAMtools/1.13-GCC-10.3.0
module load BCFtools/1.12-GCC-10.3.0

REF=/storage/research/iee_temp_dj20y461/DV_calling/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_formatted_shortnames.fna
ITER_FILE=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/sample_metadata/sample_paths_rmdup.txt
SAMPLE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $ITER_FILE | rev | cut -f1 -d'/' | rev)


WD=/storage/research/iee_temp_dj20y461/DV_calling/GLnexus/cohort_1_mt/
MT_VCF=$WD/cohort_1_chr_NC_041244.1_mitochondion_genome.vcf.gz

OUT=$WD/sample_consensus_seqs/${SAMPLE}_mt_consensus.fa

#samtools faidx $REF NC_041244.1_mitochondion_genome  |  \
#bcftools consensus -I -s $SAMPLE $MT_VCF -o /storage/research/iee_temp_dj20y461/DV_calling/GLnexus/cohort_1_mt/mt_ref.fa

#samtools faidx $REF NC_041244.1_mitochondion_genome | bcftools consensus $MT_VCF -o $OUT

samtools faidx $REF NC_041244.1_mitochondion_genome | bcftools consensus -s $SAMPLE $MT_VCF > $OUT

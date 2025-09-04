#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --export=NONE
#SBATCH --array=1-5
#SBATCH --job-name=Create_pedigree_VCFs_GATK
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training
CROSSES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/crossIDs.txt

CROSS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $CROSSES)

FATHER="${CROSS}_male_par"
MOTHER="${CROSS}_fem_par"
OFFSPRING="${CROSS}_male_1"

###################################################################
## Filter the samples for each cross from the master VCF from GATK
###############################################################

UNFILTERED_JOINT_VCF=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training/4_GATK/Unfiltered_VCF/Joint_all_unfiltered.vcf.gz
OUTDIR=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/Mendelian_evals/GATK
PEDIGREE_VCF=$OUTDIR/${CROSS}_pedigree_unfiltered.vcf.gz
MEND_EVAL_REGIONS=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/Training_genome_partitions/${CROSS}_test_partitions.bed

if [ ! -d "$OUTDIR" ]; then
   mkdir -p $OUTDIR
fi

bcftools view      $UNFILTERED_JOINT_VCF \
		-s $FATHER,$MOTHER,$OFFSPRING \
		-R $MEND_EVAL_REGIONS \
		-e 'GT[*] = "mis"' | \
bcftools view   -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' \
		-O z \
		> $PEDIGREE_VCF



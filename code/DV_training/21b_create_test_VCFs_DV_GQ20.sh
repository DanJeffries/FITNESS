#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=50G
#SBATCH --export=NONE
#SBATCH --array=1-5
#SBATCH --job-name=Create_pedigree_VCFs_DV
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

################################################################################
## Combine the individual VCFs from DV, and then filter the same way as for GATK
################################################################################

DV_VCF_FATHER=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/test/${CROSS}_male_par_test_set.vcf.gz
DV_VCF_MOTHER=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/test/${CROSS}_fem_par_test_set.vcf.gz
DV_VCF_OFFSRPING=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/test/${CROSS}_male_1_test_set.vcf.gz

JOINT_DV_VCF=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/test/${CROSS}_family.vcf.gz
JOINT_DV_VCF_FILTERED=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/test/${CROSS}_family_filtered_GQ20.vcf.gz

MEND_EVAL_REGIONS=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/Training_genome_partitions/${CROSS}_test_partitions.bed


bcftools merge $DV_VCF_FATHER $DV_VCF_MOTHER $DV_VCF_OFFSRPING -O z > $JOINT_DV_VCF

bcftools index -t $JOINT_DV_VCF

bcftools view   $JOINT_DV_VCF \
		-s $FATHER,$MOTHER,$OFFSPRING \
		-R $MEND_EVAL_REGIONS \
		-e 'GT[*] = "mis"' \
		-O z | \
bcftools view   - \
		-e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' | \
bcftools view   - \
		-i 'QUAL >= 30' \
		-O z \
		> $JOINT_DV_VCF_FILTERED



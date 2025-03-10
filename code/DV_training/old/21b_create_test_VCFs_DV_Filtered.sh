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

UNFILTERED_JOINT_VCF=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/test/${CROSS}_family.vcf.gz
MEND_EVAL_VCF_DIR=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/Mendelian_evals/DV/
FILTERED_PEDIGREE_VCF=$MEND_EVAL_VCF_DIR/${CROSS}_pedigree_ALL_FILTERS.vcf.gz
MEND_EVAL_REGIONS=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/Training_genome_partitions/${CROSS}_test_partitions.bed

if [ ! -d "$MEND_EVAL_VCF_DIR" ]; then
   mkdir $MEND_EVAL_VCF_DIR
fi

## Make filtered father VCF

FATHER_FILTERED_VCF=$MEND_EVAL_VCF_DIR/${FATHER}.vcf.gz

bcftools view      $UNFILTERED_JOINT_VCF \
                -s $FATHER \
		-R $MEND_EVAL_REGIONS \
		-O u | \
bcftools view   -i 'MIN(GQ) >= 30 & MIN(FMT/DP)>=20 & MAX(FMT/DP) <= 100' \
                -O u | \
bcftools view   -e 'GT="het" & FORMAT/AD[0:0]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.85' \
                -O u | \
bcftools view   -e 'GT="het" & FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.85' \
                -O z \
                > $FATHER_FILTERED_VCF

tabix $FATHER_FILTERED_VCF

## Make filtered mother VCF

MOTHER_FILTERED_VCF=$MEND_EVAL_VCF_DIR/${MOTHER}.vcf.gz

bcftools view      $UNFILTERED_JOINT_VCF \
                -s $MOTHER \
                -R $MEND_EVAL_REGIONS \
                -O u | \
bcftools view   -i 'MIN(GQ) >= 30 & MIN(FMT/DP)>=20 & MAX(FMT/DP) <= 100' \
                -O u | \
bcftools view   -e 'GT="het" & FORMAT/AD[0:0]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.85' \
                -O u | \
bcftools view   -e 'GT="het" & FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.85' \
                -O z \
                > $MOTHER_FILTERED_VCF

tabix $MOTHER_FILTERED_VCF


## Make filtered offspring VCF

OFFSPRING_FILTERED_VCF=$MEND_EVAL_VCF_DIR/${OFFSPRING}.vcf.gz

bcftools view      $UNFILTERED_JOINT_VCF \
                -s $OFFSPRING \
                -R $MEND_EVAL_REGIONS \
                -O u | \
bcftools view   -i 'MIN(GQ) >= 30 & MIN(FMT/DP)>=20 & MAX(FMT/DP) <= 100' \
                -O u | \
bcftools view   -e 'GT="het" & FORMAT/AD[0:0]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.85' \
                -O u | \
bcftools view   -e 'GT="het" & FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.85' \
                -O z \
                > $OFFSPRING_FILTERED_VCF

tabix $OFFSPRING_FILTERED_VCF

## Merge mother, father and offspring VCFs

bcftools merge $FATHER_FILTERED_VCF \
	       $MOTHER_FILTERED_VCF \
	       $OFFSPRING_FILTERED_VCF \
	       -O u | \
bcftools view  -e 'GT[*] = "mis"' \
               -O u | \
bcftools view  -e 'COUNT(GT="AA")=N_SAMPLES || COUNT(GT="RR")=N_SAMPLES' \
               -O z \
	       > $FILTERED_PEDIGREE_VCF

tabix $FILTERED_PEDIGREE_VCF


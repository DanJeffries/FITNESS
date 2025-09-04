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
## Concatenate the chromosome GVCF files into one GVCF per sample ##############
################################################################################

GVCFs=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training/4_GATK/GVCF

FATHER_GVCF=$GVCFs/$FATHER.all_chroms.g.vcf.gz
MOTHER_GVCF=$GVCFs/$MOTHER.all_chroms.g.vcf.gz
OFFSPRING_GVCF=$GVCFs/$OFFSPRING.all_chroms.g.vcf.gz

#bcftools concat $GVCFs/${FATHER}_NC_*vcf.gz -O z > $FATHER_GVCF
#bcftools concat $GVCFs/${MOTHER}_NC_*vcf.gz -O z > $MOTHER_GVCF
#bcftools concat $GVCFs/${OFFSPRING}_NC_*vcf.gz -O z > $OFFSPRING_GVCF

#tabix $FATHER_GVCF
#tabix $MOTHER_GVCF
#tabix $OFFSPRING_GVCF

##############################################
# make a pedigree GVCF with unfiltered loci ##
##############################################

PEDIGREE_GVCF_NO_MISSING_TEMP=$GVCFs/${CROSS}_mend_eval_chroms.NO_MISSING.TEMP.g.vcf.gz
PEDIGREE_GVCF_NO_MISSING=$GVCFs/${CROSS}_mend_eval_chroms.NO_MISSING.g.vcf.gz
MEND_EVAL_REGIONS=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/Training_genome_partitions/${CROSS}_test_partitions.bed

bcftools merge $FATHER_GVCF $MOTHER_GVCF $OFFSPRING_GVCF \
	       -O z > $PEDIGREE_GVCF_NO_MISSING_TEMP

tabix $PEDIGREE_GVCF_NO_MISSING_TEMP

bcftools view  $PEDIGREE_GVCF_NO_MISSING_TEMP \
	       -R $MEND_EVAL_REGIONS \
	       -e 'COUNT(GT="mis")>0 | COUNT(GT="RR") = 3' \
	       -a \
	       -O z > $PEDIGREE_GVCF_NO_MISSING
tabix $PEDIGREE_GVCF_NO_MISSING

#rm $PEDIGREE_GVCF_NO_MISSING_TEMP

#################################
## Filter the individual GVCFs ##
#################################

MEND_EVAL_REGIONS=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/Training_genome_partitions/${CROSS}_test_partitions.bed

FATHER_FILTERED_GVCF=$GVCFs/$FATHER.all_chroms.FILTERED.g.vcf.gz

#bcftools view      $FATHER_GVCF \
#		-i 'GT="RA" | GT="AA"' \
#                -O u | \
#bcftools view   -i 'MIN(GQ) >= 30 & MIN(FMT/DP)>=20 & MAX(FMT/DP) <= 100' \
#                -O u | \
#bcftools view   -e 'GT="het" & FORMAT/AD[0:0]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.85' \
#                -O u | \
#bcftools view   -e 'GT="het" & FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.85' \
#                -O z \
#                > $FATHER_FILTERED_GVCF

#tabix $FATHER_FILTERED_GVCF

## Make filtered mother VCF

MOTHER_FILTERED_GVCF=$GVCFs/$MOTHER.all_chroms.FILTERED.g.vcf.gz

#bcftools view     $MOTHER_GVCF \
# 		-i 'GT="RA" | GT="AA"' \
#                -O u | \
#bcftools view   -i 'MIN(GQ) >= 30 & MIN(FMT/DP)>=20 & MAX(FMT/DP) <= 100' \
#                -O u | \
#bcftools view   -e 'GT="het" & FORMAT/AD[0:0]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.85' \
#                -O u | \
#bcftools view   -e 'GT="het" & FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.85' \
#                -O z \
#               > $MOTHER_FILTERED_GVCF

#tabix $MOTHER_FILTERED_GVCF

## Make filtered offspring VCF

OFFSPRING_FILTERED_GVCF=$GVCFs/$OFFSPRING.all_chroms.FILTERED.g.vcf.gz

#bcftools view      $OFFSPRING_GVCF \
#		-i 'GT="RA" | GT="AA"' \
#                -O u | \
#bcftools view   -i 'MIN(GQ) >= 30 & MIN(FMT/DP)>=20 & MAX(FMT/DP) <= 100' \
#                -O u | \
#bcftools view   -e 'GT="het" & FORMAT/AD[0:0]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.85' \
#                -O u | \
#bcftools view   -e 'GT="het" & FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.85' \
#                -O z \
#                > $OFFSPRING_FILTERED_GVCF

#tabix $OFFSPRING_FILTERED_GVCF

###########################################
## Combine the individual Filtered GVCFs ##
###########################################

PEDIGREE_GVCF_TEMP=$GVCFs/${CROSS}_all_chroms.FILTERED_TEMP.g.vcf.gz
PEDIGREE_GVCF_FILTERED=$GVCFs/${CROSS}_all_chroms.FILTERED.MEND_EVAL.g.vcf.gz

#bcftools merge $FATHER_FILTERED_GVCF $MOTHER_FILTERED_GVCF $OFFSPRING_FILTERED_GVCF -O z > $PEDIGREE_GVCF_TEMP
#
#tabix $PEDIGREE_GVCF_TEMP
#
#bcftools view   $PEDIGREE_GVCF \
#		-R $MEND_EVAL_REGIONS \
#		-e 'GT[*] = "mis"' \
#		-O z \
#		> $PEDIGREE_GVCF_FILTERED

#tabix $PEDIGREE_GVCF_FILTERED

#rm $PEDIGREE_GVCF_TEMP

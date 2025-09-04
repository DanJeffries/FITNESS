#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --array=1
#SBATCH --job-name=Make_CONF_HOM_REF_ALT
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0
module load Anaconda3

eval "$(conda shell.bash hook)"
conda activate gvcf2bed ## for converting GVCFs to bed after filtering 

bcftools=/storage/homefs/dj20y461/Software/bcftools-1.21/bcftools


## The purpose of this script is to identify non-variant confident regions in the offspring based on the evidence in the parents. The logic is as follows:

###############################
## >>> Set paths & vars  <<< ##
###############################

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training
RESEARCH=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training
GVCFs=$RESEARCH/4_GATK/GVCF
CONF_VARS_DIR=$WD/CONF_VARS
CONFIDENT_REGIONS=$WD/Confident_regions_ALT
CROSSES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/crossIDs.txt

if [ ! -d "$CONF_VARS_DIR" ]; then
   mkdir $CONF_VARS_DIR
fi

CROSS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $CROSSES)

FATHER="${CROSS}_male_par"
MOTHER="${CROSS}_fem_par"
OFFSPRING="${CROSS}_male_1"

######################################################
## Step 1. find CONF VAR positions in the parents
######################################################

## Strategy:

#1. For parents 
#
#    - GT = DONT FILTER BY GATK GENOTYPE!??! We don't care what GATK thinks and it is super keen to call variants
#    - Min DP = 12
#    - Min ALT reads = 10 
#    - Min ALT alleles fraction = 0.8
#
#   ***NOTE - I am focusing only on Bi-allelic sites. 

MIN_DP=12
MIN_ALT_READS=10
MIN_ALT_FRACTION=0.83

## FATHER ---------------------------------------------------------------

FATHER_GVCF=$GVCFs/${FATHER}.all_chroms.g.vcf.gz
FATHER_CONF_HOM_ALT_GVCF=$GVCFs/${FATHER}.all_chroms.CONF_HOM_ALT.g.vcf.gz
FATHER_CONF_HOM_ALT_BED=$CONFIDENT_REGIONS/${FATHER}_CONF_HOM_ALT.bed

#bcftools view   $FATHER_GVCF \
#                -e 'FMT/GT="ref"' \
#                -a \
#                -O u | \
#bcftools view	-M2 \
#		-m2 \
#		-O u | \
#bcftools view   -e "FMT/DP[0] < ${MIN_DP}"  \
#		-O u | \
#bcftools view   -e "FMT/AD[0:1] < ${MIN_ALT_READS}" \
#                -O u | \
#bcftools view   -e "FMT/AD[0:1] / FMT/DP[0] < ${MIN_ALT_FRACTION}" \
#                -O z  \
#		> $FATHER_CONF_HOM_ALT_GVCF
#
#tabix $FATHER_CONF_HOM_ALT_GVCF
#
#gvcf2bed -I $FATHER_CONF_HOM_ALT_GVCF -O $FATHER_CONF_HOM_ALT_BED -q 0 -nq 0

## MOTHER ---------------------------------------------------------------

MOTHER_GVCF=$GVCFs/${MOTHER}.all_chroms.g.vcf.gz
MOTHER_CONF_HOM_ALT_GVCF=$GVCFs/${MOTHER}.all_chroms.CONF_HOM_ALT.g.vcf.gz
MOTHER_CONF_HOM_ALT_BED=$CONFIDENT_REGIONS/${MOTHER}_CONF_HOM_ALT.bed

#bcftools view   $MOTHER_GVCF \
#                -e 'FMT/GT="ref"' \
#                -a \
#                -O u | \
#bcftools view   -M2 \
#                -m2 \
#                -O u | \
#bcftools view   -e "FMT/DP[0] < ${MIN_DP}"  \
#                -O u | \
#bcftools view   -e "FMT/AD[0:1] < ${MIN_ALT_READS}" \
#                -O u | \
#bcftools view   -e "FMT/AD[0:1] / FMT/DP[0] < ${MIN_ALT_FRACTION}" \
#                -O z  \
#                > $MOTHER_CONF_HOM_ALT_GVCF
#
#tabix $MOTHER_CONF_HOM_ALT_GVCF
#
#gvcf2bed -I $MOTHER_CONF_HOM_ALT_GVCF -O $MOTHER_CONF_HOM_ALT_BED -q 0 -nq 0

## OFFSPRING ---------------------------------------------------------------

OFFSPRING_GVCF=$GVCFs/${OFFSPRING}.all_chroms.g.vcf.gz
OFFSPRING_CONF_GVCF=$GVCFs/${OFFSPRING}.all_chroms.CONF_ALT.g.vcf.gz
OFFSPRING_CONF_BED=$CONFIDENT_REGIONS/${OFFSPRING}_conf_ALT.bed

###########################################################################
# Finding CONF VARs in offspring ------------

FATHER_CONF_HOM_REF_BED=$CONFIDENT_REGIONS/${FATHER}_conf_ALT.bed
MOTHER_CONF_HOM_REF_BED=$CONFIDENT_REGIONS/${MOTHER}_conf_ALT.bed
OFFSPRING_HET_OO11_BED=${CONF_VARS_DIR}/${OFFSPRING}_CONF_HET_0011.bed
OFFSPRING_HET_1100_BED=${CONF_VARS_DIR}/${OFFSPRING}_CONF_HET_1100.bed
OFFSPRING_HOM_ALT_BED=${CONF_VARS_DIR}/${OFFSPRING}_CONF_HOM_ALT.bed

## 0/0 father and 1/1 mother ------------------------------------------------

OFFSPRING_HET_0011_VCF=$GVCFs/${OFFSPRING}.all_chroms.HET_0011.g.vcf.gz

#bedtools intersect -a $FATHER_CONF_HOM_REF_BED -b $MOTHER_CONF_HOM_ALT_BED > $OFFSPRING_HET_OO11_BED

MIN_DP=10
MIN_ALLELE_FRACTION=0.2

#bcftools view     $OFFSPRING_GVCF \
#		-R $OFFSPRING_HET_OO11_BED \
#		-a \
#		-i "FMT/GT='het'" \
#		-O u | \
#bcftools view   -e "FMT/DP<$MIN_DP | FMT/AD[0:0] < ${MIN_ALLELE_FRACTION} | FMT/AD[0:1] < ${MIN_ALLELE_FRACTION}" \
#                -O z \
#		>  $OFFSPRING_HET_0011_VCF
#
#tabix $OFFSPRING_HET_0011_VCF

## 1/1 father and 0/0 mother ------------------------------------------------

OFFSPRING_HET_1100_VCF=$GVCFs/${OFFSPRING}.all_chroms.HET_1100.g.vcf.gz

#bedtools intersect -a $FATHER_CONF_HOM_ALT_BED -b $MOTHER_CONF_HOM_REF_BED > $OFFSPRING_HET_1100_BED

#bcftools view     $OFFSPRING_GVCF \
#		-R $OFFSPRING_HET_1100_BED \
#		-a \
#		-i "FMT/GT='het'" \
#		-O u | \
#bcftools view  -e "FMT/DP<$MIN_DP | FMT/AD[0:0] < ${MIN_ALLELE_FRACTION} | FMT/AD[0:1] < ${MIN_ALLELE_FRACTION}" \
#		-O z \
#		>  $OFFSPRING_HET_1100_VCF
#
#tabix $OFFSPRING_HET_1100_VCF
#
## 1/1 father and 1/1 mother -------------------------------------------------

OFFSPRING_HOM_ALT_VCF=$GVCFs/${OFFSPRING}.all_chroms.HOM_ALT.g.vcf.gz

#bedtools intersect -a $FATHER_CONF_HOM_ALT_BED -b $MOTHER_CONF_HOM_ALT_BED > $OFFSPRING_HOM_ALT_BED

bcftools view     $OFFSPRING_GVCF \
		-R $OFFSPRING_HOM_ALT_BED \
		-a \
		-i "FMT/GT='AA'" \
		-O u | \
bcftools view  -e "FMT/DP<$MIN_DP | FMT/AD[0:1] / FMT/DP[0] < ${MIN_ALT_FRACTION}" \
		-O z \
		>  $OFFSPRING_HOM_ALT_VCF

tabix $OFFSPRING_HOM_ALT_VCF

### Finally, merge the VCFs and convert to BED

OFFSPRING_CONF_VARS_VCF=$GVCFs/${OFFSPRING}.all_chroms.CONF_VARS.vcf.gz
OFFSPRING_CONF_VARS_BED=${CONF_VARS_DIR}/${OFFSPRING}_CONF_VARS_ALL.bed

bcftools concat $OFFSPRING_HET_0011_VCF \
	        $OFFSPRING_HET_1100_VCF \
		$OFFSPRING_HOM_ALT_VCF \
		-a \
		-O z \
		> $OFFSPRING_CONF_VARS_VCF

tabix $OFFSPRING_CONF_VARS_VCF

gvcf2bed -I $OFFSPRING_CONF_VARS_VCF -O $OFFSPRING_CONF_VARS_BED -q 0 -nq 0 


############################################################################
## NEXT ADD THE CONF VAR REGIONS TO THE HOM_REF REGIONS
#######################################################





######################################################
## Step 3. make the mask
######################################################

GA_1N_WINDOW_MASK=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training/5_masks/Finding_2n_windows/Ga_1n_windows.bed ## 1n windows on the sex chromosome
REPEAT_MASK=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training/5_masks/repeats_renamed.bed ## custom rep library
COMBINED_MASK=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training/5_masks/combined_mask.bed 

## combine the repeats and 1n mask

#cat $GA_1N_WINDOW_MASK $REPEAT_MASK | bedtools sort | bedtools merge > $COMBINED_MASK  ## combine repeat and 1n window mask

###################################################################
## Step 4. Remove the masked regions from the CONF_REGIONS BED file
###################################################################

CONF_REGIONS_MASKED_BED=$CONFIDENT_REGIONS/${CROSS}_CONF_HOM_REF_ALT_masked.bed

#bedtools subtract -a $OFFSPRING_CONF_BED_SHARED -b $COMBINED_MASK | \
#bedtools sort | \
#bedtools merge 	> $CONF_REGIONS_MASKED_BED



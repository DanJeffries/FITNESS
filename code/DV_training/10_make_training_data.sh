#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=08:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --array=1-5
#SBATCH --job-name=Make_training_data
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0
module load Anaconda3

eval "$(conda shell.bash hook)"
conda activate gvcf2bed ## for converting GVCFs to bed after filtering 

bcftools=/storage/homefs/dj20y461/Software/bcftools-1.21/bcftools

###############################
## >>> Set paths & variables  <<< ##
###############################

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/TRAINING_DATA
GVCFs=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training/4_GATK/GVCF
CROSSES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/crossIDs.txt

CROSS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $CROSSES)

FATHER="${CROSS}_male_par"
MOTHER="${CROSS}_fem_par"
OFFSPRING="${CROSS}_male_1"

FATHER_GVCF=$GVCFs/${FATHER}.all_chroms.g.vcf.gz
MOTHER_GVCF=$GVCFs/${MOTHER}.all_chroms.g.vcf.gz
OFFSPRING_GVCF=$GVCFs/${OFFSPRING}.all_chroms.g.vcf.gz

############################################################################################################
############################  FIND CONFIDENT HOM_REF REGIONS  ##############################################
############################################################################################################

######################################################
## Step 1. find 0/0 confident sites in each sample
######################################################

### CONF HOM_REF

#    - GT = DONT FILTER BY GATK GENOTYPE!??! We don't care what GATK thinks and it is super keen to call variants. E.g. 17,3 will be called 0/1 in GATK.
#    - Min DP = 15
#    - Min ref reads = 13
#    - Min REF alleles fraction = 0.86
#

echo "Making CONF HOM_REF" 

MIN_DP=15
MIN_REF_READS=13
MIN_REF_FRACTION=0.86

## FATHER ---------------------------------------------------------------

FATHER_CONF_HOM_REF_GVCF=$WD/${FATHER}.CONF_HOM_REF.g.vcf.gz
FATHER_CONF_HOM_REF_BED=$WD/${FATHER}.CONF_HOM_REF.bed

bcftools view   $FATHER_GVCF \
                -e "FMT/DP[0] < ${MIN_DP}" \
                -a \
                -O u | \
bcftools view   -e "FMT/AD[0:0] < ${MIN_REF_READS}" \
                -O u | \
bcftools view   -e "(FMT/AD[0:0] / FMT/DP[0]) < ${MIN_REF_FRACTION}" \
                -O z \
                > $FATHER_CONF_HOM_REF_GVCF

tabix $FATHER_CONF_HOM_REF_GVCF

gvcf2bed -I $FATHER_CONF_HOM_REF_GVCF -O $FATHER_CONF_HOM_REF_BED -q 0 -nq 0

echo "   FATHER CONF HOM REF FILES:"
echo "          $FATHER_CONF_HOM_REF_GVCF"
echo "          $FATHER_CONF_HOM_REF_BED"

## MOTHER ---------------------------------------------------------------

MOTHER_CONF_HOM_REF_GVCF=$WD/${MOTHER}.CONF_HOM_REF.g.vcf.gz
MOTHER_CONF_HOM_REF_BED=$WD/${MOTHER}.CONF_HOM_REF.bed

bcftools view   $MOTHER_GVCF \
                -e "FMT/DP[0] < ${MIN_DP}" \
                -a \
                -O u | \
bcftools view   -e "FMT/AD[0:0] < ${MIN_REF_READS}" \
                -O u | \
bcftools view   -e "(FMT/AD[0:0] / FMT/DP[0]) < ${MIN_REF_FRACTION}" \
		-O z \
		> $MOTHER_CONF_HOM_REF_GVCF

tabix $MOTHER_CONF_HOM_REF_GVCF

gvcf2bed -I $MOTHER_CONF_HOM_REF_GVCF -O $MOTHER_CONF_HOM_REF_BED -q 0 -nq 0

echo "   MOTHER CONF HOM REF FILES:"
echo "          $MOTHER_CONF_HOM_REF_GVCF"
echo "          $MOTHER_CONF_HOM_REF_BED"

## OFFSPRING ---------------------------------------------------------------

OFFSPRING_FILTERED_HOM_REF_GVCF=$WD/${OFFSPRING}.FILTERED_HOM_REF.g.vcf.gz
OFFSPRING_FILTERED_HOM_REF_BED=$WD/${OFFSPRING}.FILTERED_HOM_REF.bed

bcftools view   $OFFSPRING_GVCF \
                -e "FMT/DP[0] < ${MIN_DP}" \
                -a \
                -O u | \
bcftools view   -e "FMT/AD[0:0] < ${MIN_REF_READS}" \
                -O u | \
bcftools view   -e "(FMT/AD[0:0] / FMT/DP[0]) < ${MIN_REF_FRACTION}" \
                -O z \
                > $OFFSPRING_FILTERED_HOM_REF_GVCF

tabix $OFFSPRING_FILTERED_HOM_REF_GVCF

gvcf2bed -I $OFFSPRING_FILTERED_HOM_REF_GVCF -O $OFFSPRING_FILTERED_HOM_REF_BED -q 0 -nq 0

echo "   OFFSPRING CONF HOM REF FILES:"
echo "          $OFFSPRING_FILTERED_HOM_REF_GVCF"
echo "          $OFFSPRING_FILTERED_HOM_REF_BED"

#################################################################################
## Step 2. Keep only CONF HOM_REF positions found in offspring, mother and father.
#################################################################################

PARENTS_SHARED_CONF_HOM_REF_BED=$WD/${CROSS}_parents.SHARED_CONF_HOM_REF.bed
OFFSPRING_CONF_HOM_REF_BED=$WD/${OFFSPRING}.CONF_HOM_REF.bed

bedtools intersect -a $FATHER_CONF_HOM_REF_BED -b $MOTHER_CONF_HOM_REF_BED  > $PARENTS_SHARED_CONF_HOM_REF_BED
bedtools intersect -a $OFFSPRING_FILTERED_HOM_REF_BED -b $PARENTS_SHARED_CONF_HOM_REF_BED > $OFFSPRING_CONF_HOM_REF_BED

###################################################################################
################# 	MAKE CONF VARIANTS 	###################################
###################################################################################

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

FATHER_CONF_HOM_ALT_GVCF=$WD/${FATHER}.CONF_HOM_ALT.g.vcf.gz
FATHER_CONF_HOM_ALT_BED=$WD/${FATHER}.CONF_HOM_ALT.bed

bcftools view   $FATHER_GVCF \
                -e 'FMT/GT="ref"' \
                -a \
                -O u | \
bcftools view	-M2 \
		-m2 \
		-O u | \
bcftools view   -e "FMT/DP[0] < ${MIN_DP}"  \
		-O u | \
bcftools view   -e "FMT/AD[0:1] < ${MIN_ALT_READS}" \
                -O u | \
bcftools view   -e "FMT/AD[0:1] / FMT/DP[0] < ${MIN_ALT_FRACTION}" \
                -O z  \
		> $FATHER_CONF_HOM_ALT_GVCF

tabix $FATHER_CONF_HOM_ALT_GVCF
gvcf2bed -I $FATHER_CONF_HOM_ALT_GVCF -O $FATHER_CONF_HOM_ALT_BED -q 0 -nq 0

echo "   FATHER CONF HOM ALT FILES:"
echo "          $FATHER_CONF_HOM_ALT_GVCF"
echo "          $FATHER_CONF_HOM_ALT_BED"

## MOTHER ---------------------------------------------------------------

MOTHER_CONF_HOM_ALT_GVCF=$WD/${MOTHER}.CONF_HOM_ALT.g.vcf.gz
MOTHER_CONF_HOM_ALT_BED=$WD/${MOTHER}.CONF_HOM_ALT.bed

bcftools view   $MOTHER_GVCF \
                -e 'FMT/GT="ref"' \
                -a \
                -O u | \
bcftools view   -M2 \
                -m2 \
                -O u | \
bcftools view   -e "FMT/DP[0] < ${MIN_DP}"  \
                -O u | \
bcftools view   -e "FMT/AD[0:1] < ${MIN_ALT_READS}" \
                -O u | \
bcftools view   -e "FMT/AD[0:1] / FMT/DP[0] < ${MIN_ALT_FRACTION}" \
                -O z  \
                > $MOTHER_CONF_HOM_ALT_GVCF

tabix $MOTHER_CONF_HOM_ALT_GVCF
gvcf2bed -I $MOTHER_CONF_HOM_ALT_GVCF -O $MOTHER_CONF_HOM_ALT_BED -q 0 -nq 0

echo "   MOTHER CONF HOM ALT FILES:"
echo "          $MOTHER_CONF_HOM_ALT_GVCF"
echo "          $MOTHER_CONF_HOM_ALT_BED"

## OFFSPRING ---------------------------------------------------------------

## 0/0 father and 1/1 mother ------------------------------------------------

OFFSPRING_HET_0011_VCF=$WD/${OFFSPRING}.CONF_HET_0011.g.vcf.gz
OFFSPRING_HET_OO11_BED=$WD/${OFFSPRING}.CONF_HET_0011.bed

bedtools intersect -a $FATHER_CONF_HOM_REF_BED -b $MOTHER_CONF_HOM_ALT_BED > $OFFSPRING_HET_OO11_BED

MIN_DP=10
MIN_ALLELE_FRACTION=0.2

bcftools view     $OFFSPRING_GVCF \
		-R $OFFSPRING_HET_OO11_BED \
		-a \
		-i "FMT/GT='het'" \
		-O u | \
bcftools view   -e "FMT/DP<$MIN_DP | FMT/AD[0:0] < ${MIN_ALLELE_FRACTION} | FMT/AD[0:1] < ${MIN_ALLELE_FRACTION}" \
                -O z \
		>  $OFFSPRING_HET_0011_VCF

tabix $OFFSPRING_HET_0011_VCF

## 1/1 father and 0/0 mother ------------------------------------------------

OFFSPRING_HET_1100_VCF=$WD/${OFFSPRING}.CONF_HET_1100.g.vcf.gz
OFFSPRING_HET_1100_BED=$WD/${OFFSPRING}.CONF_HET_1100.bed

bedtools intersect -a $FATHER_CONF_HOM_ALT_BED -b $MOTHER_CONF_HOM_REF_BED > $OFFSPRING_HET_1100_BED

bcftools view     $OFFSPRING_GVCF \
		-R $OFFSPRING_HET_1100_BED \
		-a \
		-i "FMT/GT='het'" \
		-O u | \
bcftools view  -e "FMT/DP<$MIN_DP | FMT/AD[0:0] < ${MIN_ALLELE_FRACTION} | FMT/AD[0:1] < ${MIN_ALLELE_FRACTION}" \
		-O z \
		>  $OFFSPRING_HET_1100_VCF

tabix $OFFSPRING_HET_1100_VCF
#
## 1/1 father and 1/1 mother -------------------------------------------------

OFFSPRING_HOM_ALT_VCF=$WD/${OFFSPRING}.CONF_HOM_ALT.g.vcf.gz
OFFSPRING_HOM_ALT_BED=$WD/${OFFSPRING}.CONF_HOM_ALT.bed

bedtools intersect -a $FATHER_CONF_HOM_ALT_BED -b $MOTHER_CONF_HOM_ALT_BED > $OFFSPRING_HOM_ALT_BED

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

OFFSPRING_CONF_VARS_VCF=$WD/${OFFSPRING}.CONF_VARS_ALL.vcf.gz
OFFSPRING_CONF_VARS_BED=$WD/${OFFSPRING}.CONF_VARS_ALL.bed

bcftools concat $OFFSPRING_HET_0011_VCF \
	        $OFFSPRING_HET_1100_VCF \
		$OFFSPRING_HOM_ALT_VCF \
		-a \
		-O z \
		> $OFFSPRING_CONF_VARS_VCF

tabix $OFFSPRING_CONF_VARS_VCF

gvcf2bed -I $OFFSPRING_CONF_VARS_VCF -O $OFFSPRING_CONF_VARS_BED -q 0 -nq 0 

echo "   OFFSPRING CONFIDENT VARIANT FILES:"
echo "          $OFFSPRING_CONF_VARS_VCF"
echo "          $OFFSPRING_CONF_VARS_BED"


############################################################################
## NEXT ADD THE CONF VAR REGIONS TO THE HOM_REF REGIONS
#######################################################

OFFSPRING_CONF_REGIONS_BED=$WD/${OFFSPRING}.CONF_REGIONS.bed

cat $OFFSPRING_CONF_HOM_REF_BED $OFFSPRING_CONF_VARS_BED | bedtools sort | bedtools merge > $OFFSPRING_CONF_REGIONS_BED

###########################################################
## Step 3. remove masked regions from confident regions
###########################################################

GA_1N_WINDOW_MASK=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training/5_masks/Finding_2n_windows/Ga_1n_windows.bed ## 1n windows on the sex chromosome
REPEAT_MASK=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training/5_masks/repeats_renamed.bed ## custom rep library
COMBINED_MASK=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training/5_masks/combined_mask.bed 

## combine the repeats and 1n mask

cat $GA_1N_WINDOW_MASK $REPEAT_MASK | bedtools sort | bedtools merge > $COMBINED_MASK  ## combine repeat and 1n window mask

###################################################################
## Step 4. Remove the masked regions from the CONF_REGIONS BED file
###################################################################

OFFSPRING_CONF_REGIONS_MASKED_BED=$WD/${OFFSPRING}.CONF_REGIONS_MASKED.bed

bedtools subtract -a $OFFSPRING_CONF_REGIONS_BED -b $COMBINED_MASK | \
bedtools sort | \
bedtools merge 	> $OFFSPRING_CONF_REGIONS_MASKED_BED

echo "   OFFSPRING MASKED CONFIDENT REGIONS FILE HERE:"
echo "          $OFFSPRING_CONF_REGIONS_MASKED_BED"

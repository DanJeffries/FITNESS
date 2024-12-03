#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --array=1-5
#SBATCH --job-name=Make_non_var_conf_regions
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0

bcftools=/storage/homefs/dj20y461/Software/bcftools-1.21/bcftools

## The purpose of this script is to identify non-variant confident regions in the offspring based on the evidence in the parents. The logic is as follows:

# In order to train the model to learn what a non-variant position looks like, the "confident regions" must contain loci which we are confident are non-variant in the offspring. However, we do not want to use
# evidence in the offspring to identify these positions. For example, if we removed loci called as 0/1 in the offspring (despite parents both being 0/0) we would be removing false positive variant signal. We want DV to learn what false positives
# look like, so we want these situations in there. 

# The solution is to identify genomic regions that we EXPECT to be homozygous in the offspring based on the parents genotypes. Such loci should have identical genotypes in the parents, which are called with extremely good confidence. 
# Note that it is not enough just that both parents be homozygous, they must be homozygous for the same allele. E.g. AA & AA, not AA & TT. 

# The below filters will identify such regions in the parents:

###############################
## >>> Set paths & vars  <<< ##
###############################

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training
RESEARCH=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training
GVCFs=$RESEARCH/4_GATK/GVCF
MERGED_GVCFs=$WD/GVCFs_merged
CONFIDENT_REGIONS=$WD/Confident_regions
CROSSES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/crossIDs.txt
FILTERED_GVCF_DIR=$WD/Filtered_GVCFs

if [ ! -d "$MERGED_GVCFs" ]; then
   mkdir $MERGED_GVCFs
fi

if [ ! -d "$FILTERED_GVCF_DIR" ]; then
   mkdir $FILTERED_GVCF_DIR
fi

CROSS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $CROSSES)

FATHER="${CROSS}_male_par"
MOTHER="${CROSS}_fem_par"
OFFSPRING="${CROSS}_male_1"

######################################################
## >>> CONCAT chromosome GVCFs for each samaple <<< ##
######################################################

# GVCFs are output for each chromosome from my parallelised GATK script. So we will concatenate them first

FATHER_MERGED_GVCF=$MERGED_GVCFs/$FATHER.g.vcf.gz
MOTHER_MERGED_GVCF=$MERGED_GVCFs/$MOTHER.g.vcf.gz
OFFSPRING_MERGED_GVCF=$MERGED_GVCFs/$OFFSPRING.g.vcf.gz

#bcftools concat -O z $GVCFs/$FATHER*chromosome*gz > $FATHER_MERGED_GVCF
#bcftools concat -O z $GVCFs/$MOTHER*chromosome*gz > $MOTHER_MERGED_GVCF
#bcftools concat -O z $GVCFs/$OFFSPRING*chromosome*gz > $OFFSPRING_MERGED_GVCF

#tabix $FATHER_MERGED_GVCF
#tabix $MOTHER_MERGED_GVCF
#tabix $OFFSPRING_MERGED_GVCF


########################################
### >>> Filter parent GVCFs  <<< ####
########################################

# So here I first combine the individual concatenated GVCFs of each parent into a single GVCF for each CROSS, allowing me to do some sample-wise filtering below. This will include ensuring that all genotypes are homozygous reference, minimum depth across the two parents is 20x, and the genotype quality for both parents is good GQ>=30.

PARENTS_GVCF=$MERGED_GVCFs/${CROSS}_parents.g.vcf.gz

bcftools merge $MERGED_GVCFs/${FATHER}.g.vcf.gz \
	       $MERGED_GVCFs/${MOTHER}.g.vcf.gz \
	       -O u | \
bcftools view  -i 'COUNT(GT="RR")=2 & MIN(FMT/DP)>=20 & MIN(GQ)>=30 & N_MISSING=0' \
	       -O z  \
               > $PARENTS_GVCF

tabix $PARENTS_GVCF

########################################
## >>> converting VCF coordinates <<< ##
########################################

# Then I need to convert these regions to bed format to give to DV. I will also use this to creat the offspring conf regions VCF to check what these loci look like in the offspring. 

## NOTE: VCF and BED have different indexing systems. So when we convert formats we need to convert coordinates too.

# - VCF = 1-based, closed [start, end]
# - BED = 0-based, half-open [start-1, start)

# 1-based means that the number 1 refers to the 1st element, number 2 refers to the 2nd etc
# 0-based means that 0 refers to the first element, 1 refers to the 2nd element etc.

# Closed indexing [denoted with square brackets] refers to indexes that INCLUDE the end point. So [2-4] includes the 2nd, 3rd, 4th elements (if using 1-based indexing)
# Open indexing (denoted with curved brackets) refers to indexes that EXCLUDE the end points. So (2-4) includes only the 3rd element (if using 1-based indexing)

# So to convert from VCF -> BED coordinates, we need to
# 1. convert from 1-based to 0-based indexing.
#    1. To do this we just -1 from every number
# 2. Convert from closed [] to half-open [)
#    1. So we don’t need to do anything to the start point.
#    2. But the second endpoint would not be included, so we need to +1 to this.

# So we need to -1 from both start and end point, and +1 to the end point. The latter cancels out so we do nothing to the end point.

#####################################################################################
# >> In summary, to convert VCF coords -> BED coords, we -1 from the start point. << #
#####################################################################################

CONF_REGIONS_BED=$CONFIDENT_REGIONS/${CROSS}_nonVar_conf_regions.bed

echo "starting bed creation"
zcat $PARENTS_GVCF | cut -f1,2,8  | grep 'END' | grep -v '#' | sed 's/END=//g' | cut -f1 -d';' | \
awk '{print $1 "\t" ($2 - 1) "\t" $3}' | grep -v 'BaseQRankSum' | grep -v 'ExcessHet' | bedtools merge > $CONF_REGIONS_BED

echo "bed created"

### do i need to add a missing values filter in? I think not, DV takes the bam file, so if there is no call in the offspring, it means that there were not enough reads for GATK, but there might be for DV. Eitherway, it is a real world representation of what homozygous regions can look like.  


#################################################################################
 #### FINALLY REMOVE THE REPEAT MASK AND 1N WINDOWS FROM THE CONF REGIONS ######
#################################################################################

GA_1N_WINDOW_MASK=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training/5_masks/Finding_2n_windows/Ga_1n_windows.bed
REPEAT_MASK=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training/5_masks/repeats_renamed.bed
COMBINED_MASK=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training/5_masks/combined_mask.bed

## combine the repeats and 1n mask - these are regions we want to exclude. Need to remove these from the confident regions

CROSS_CONF_BED_FINAL=$CONFIDENT_REGIONS/${CROSS}_nonVar_conf_regions_masked.bed

cat $GA_1N_WINDOW_MASK $REPEAT_MASK | bedtools sort | bedtools merge > $COMBINED_MASK  ## combine repeat and 1n window mask

## REMOVE THESE REGIONS FROM THE CONFIDENT REGIONS IDENTIFIED FROM THE PARENTS ##

bedtools subtract -a $CONF_REGIONS_BED -b $COMBINED_MASK > $CROSS_CONF_BED_FINAL  ## subtract combined repeat and 1n window mask from the confident regions

#################################################################################
################ CREATE OFFSPRING CONFIDENT_REGIONS GVCF  ### ###################
#################################################################################

# This will be useful for validation and summary plots

CROSS_CONF_REGIONS_GVCF=$CONFIDENT_REGIONS/${CROSS}_male_1_putative_nonVar_confident_regions.g.vcf.gz

$bcftools view $OFFSPRING_MERGED_GVCF \
              -R $CROSS_CONF_BED_FINAL \
	      -i "N_MISSING=0" \
              -a \
              -O z \
              > $CROSS_CONF_REGIONS_GVCF

tabix $CROSS_CONF_REGIONS_GVCF

### do i need to add a missing values filter in? I think not, DV takes the bam file, so if there is no call in the offspring, it means that there were not enough reads for GATK, but there might be for DV. Eitherway, it is a real world representation of what homozygous regions can look like.





#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --array=2-5
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

# In order to train the model to learn what a non-variant position looks like, the "confident regions" must contain loci which we are confident are non-variant in the offspring. However, we do not want to use
# evidence in the offspring to identify these positions. For example, if we removed loci called as 0/1 in the offspring (despite parents both being 0/0) we would be removing false positive variant signal. We want DV to learn what false positives
# look like, so we want these situations in there. 

# The solution is to identify genomic regions that we EXPECT to be homozygous in the offspring based on the parents genotypes. Such loci should have identical genotypes in the parents, which are called with extremely good confidence. 
# Note that it is not enough just that both parents be homozygous, they must be homozygous for the same allele. E.g. AA & AA, not AA & TT. 

########################################
##### >>> Converting from VCF to BED
########################################

## Note - converting VCF to bed coordinates is a pain . . . see below:

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

# There are extra considerations for INDELS. I.e. where the allele length is greater than 1. The question being do we add the allele length to the end position of the bed

## In the end I found a tool (gvcf2bed) that will do the conversion. It takes care of the INDEL length thing. 

###############################
## >>> Set paths & vars  <<< ##
###############################

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training
RESEARCH=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training
GVCFs=$RESEARCH/4_GATK/GVCF
CONFIDENT_REGIONS=$WD/Confident_regions_ALT
CROSSES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/crossIDs.txt

if [ ! -d "$CONFIDENT_REGIONS" ]; then
   mkdir $CONFIDENT_REGIONS
fi

CROSS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $CROSSES)

FATHER="${CROSS}_male_par"
MOTHER="${CROSS}_fem_par"
OFFSPRING="${CROSS}_male_1"

######################################################
## Step 1. find 0/0 confident sites in each parent
######################################################

## New filtering strategy. 

### CONF HOM_REF

#1. For offspring ONLY 
#
#    - GT = DONT FILTER BY GATK GENOTYPE!??! We don't care what GATK thinks and it is super keen to call variants. E.g. 17,3 will be called 0/1 in GATK.
#    - Min DP = 15
#    - Min ref reads = 13 
#    - Min REF alleles fraction = 0.8
#

MIN_DP=15
MIN_REF_READS=13
MIN_REF_FRACTION=0.86


## FATHER ---------------------------------------------------------------

FATHER_GVCF=$GVCFs/${FATHER}.all_chroms.g.vcf.gz
FATHER_CONF_GVCF=$GVCFs/${FATHER}.all_chroms.CONF_ALT.g.vcf.gz
FATHER_CONF_BED=$CONFIDENT_REGIONS/${FATHER}_conf_ALT.bed

bcftools view   $FATHER_GVCF \
                -e "FMT/DP[0] < ${MIN_DP}" \
                -a \
                -O u | \
bcftools view   -e "FMT/AD[0:0] < ${MIN_REF_READS}" \
                -O u | \
bcftools view   -e "(FMT/AD[0:0] / FMT/DP[0]) < ${MIN_REF_FRACTION}" \
                -O z \
                > $FATHER_CONF_GVCF

tabix $FATHER_CONF_GVCF

gvcf2bed -I $FATHER_CONF_GVCF -O $FATHER_CONF_BED -q 0 -nq 0

## MOTHER ---------------------------------------------------------------

MOTHER_GVCF=$GVCFs/${MOTHER}.all_chroms.g.vcf.gz
MOTHER_CONF_GVCF=$GVCFs/${MOTHER}.all_chroms.CONF_ALT.g.vcf.gz
MOTHER_CONF_BED=$CONFIDENT_REGIONS/${MOTHER}_conf_ALT.bed

bcftools view   $MOTHER_GVCF \
                -e "FMT/DP[0] < ${MIN_DP}" \
                -a \
                -O u | \
bcftools view   -e "FMT/AD[0:0] < ${MIN_REF_READS}" \
                -O u | \
bcftools view   -e "(FMT/AD[0:0] / FMT/DP[0]) < ${MIN_REF_FRACTION}" \
                -O z \
                > $MOTHER_CONF_GVCF

tabix $MOTHER_CONF_GVCF

gvcf2bed -I $MOTHER_CONF_GVCF -O $MOTHER_CONF_BED -q 0 -nq 0


## OFFSPRING ---------------------------------------------------------------

OFFSPRING_GVCF=$GVCFs/${OFFSPRING}.all_chroms.g.vcf.gz
OFFSPRING_CONF_GVCF=$GVCFs/${OFFSPRING}.all_chroms.CONF_ALT.g.vcf.gz
OFFSPRING_CONF_BED=$CONFIDENT_REGIONS/${OFFSPRING}_conf_ALT.bed

bcftools view   $OFFSPRING_GVCF \
                -e "FMT/DP[0] < ${MIN_DP}" \
                -a \
                -O u | \
bcftools view	-e "FMT/AD[0:0] < ${MIN_REF_READS}" \
		-O u | \
bcftools view	-e "(FMT/AD[0:0] / FMT/DP[0]) < ${MIN_REF_FRACTION}" \
		-O z \
                > $OFFSPRING_CONF_GVCF

tabix $OFFSPRING_CONF_GVCF

gvcf2bed -I $OFFSPRING_CONF_GVCF -O $OFFSPRING_CONF_BED -q 0 -nq 0


########################################################################
## Step 2. filter for positions kept in all samples I.e. the intersection 
#######################################################################

PARENTS_CONF_BED=$CONFIDENT_REGIONS/${CROSS}_parents_conf.bed
OFFSPRING_CONF_BED_SHARED=$CONFIDENT_REGIONS/${OFFSPRING}_conf_ALT_SHARED.bed

bedtools intersect -a $FATHER_CONF_BED -b $MOTHER_CONF_BED  > $PARENTS_CONF_BED
bedtools intersect -a $OFFSPRING_CONF_BED -b $PARENTS_CONF_BED > $OFFSPRING_CONF_BED_SHARED

######################################################
## Step 3. make the mask
######################################################

GA_1N_WINDOW_MASK=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training/5_masks/Finding_2n_windows/Ga_1n_windows.bed ## 1n windows on the sex chromosome
REPEAT_MASK=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training/5_masks/repeats_renamed.bed ## custom rep library
COMBINED_MASK=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training/5_masks/combined_mask.bed 

## combine the repeats and 1n mask

#cat $GA_1N_WINDOW_MASK $REPEAT_MASK | bedtools sort | bedtools merge > $COMBINED_MASK  ## combine repeat and 1n window mask

###################################################################
## Step 4. Remove the masked regions from the confident bed file
###################################################################

CONF_REGIONS_MASKED_BED=$CONFIDENT_REGIONS/${CROSS}_CONF_HOM_REF_ALT_masked.bed

bedtools subtract -a $OFFSPRING_CONF_BED_SHARED -b $COMBINED_MASK | \
bedtools sort | \
bedtools merge 	> $CONF_REGIONS_MASKED_BED



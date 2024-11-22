#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --export=NONE
#SBATCH --array=1-5
#SBATCH --job-name=GVCF_2_BED
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0
module load BEDTools/2.30.0-GCC-10.3.0

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training
CROSSES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/crossIDs.txt
FILTERED_GVCF_DIR=$WD/Filtered_GVCFs

CONF_BED_DIR=$WD/Confident_regions

if [ ! -d "$CONF_BED_DIR" ]; then
   mkdir $CONF_BED_DIR
fi

CROSS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $CROSSES)

FATHER_CONF_GVCF=$FILTERED_GVCF_DIR/${CROSS}_male_par.conf.g.vcf.gz
MOTHER_CONF_GVCF=$FILTERED_GVCF_DIR/${CROSS}_fem_par.conf.g.vcf.gz

FATHER_BED=$CONF_BED_DIR/${CROSS}_male_par.conf.bed
MOTHER_BED=$CONF_BED_DIR/${CROSS}_fem_par.conf.bed

OFFSPRING_BED=$CONF_BED_DIR/${CROSS}_male_1.conf.bed

##################################################
## converting VCF coordinates
############################

## VCF and BED have different indexing systems. So when we convert formats we need to convert coordinates too. 

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


## FATHER
zcat $FATHER_CONF_GVCF | cut -f1,2,8  | grep 'END' | grep -v '#' | sed 's/END=//g' | \
awk '{print $1 "\t" ($2 - 1) "\t" $3}' | bedtools merge > $FATHER_BED

## MOTHER
zcat $MOTHER_CONF_GVCF | cut -f1,2,8  | grep 'END' | grep -v '#' | sed 's/END=//g' | \
awk '{print $1 "\t" ($2 - 1) "\t" $3}' | bedtools merge > $MOTHER_BED

#################################################################################
#### REMOVE THE REPEAT MASK AND 1N WINDOWS FROM THE CONF REGIONS ######
#################################################################################

GA_1N_WINDOW_MASK=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training/Finding_2n_windows/Ga_1n_windows.bed
REPEAT_MASK=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training/masks/repeats_renamed.bed
COMBINED_MASK=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training/masks/combined_mask.bed

## combine the repeats and 1n mask - these are regions we want to exclude. Need to remove these from the confident regions 

cat $GA_1N_WINDOW_MASK $REPEAT_MASK | bedtools sort | bedtools merge > $COMBINED_MASK

## REMOVE THESE REGIONS FROM THE CONFIDENT REGIONS IN FATHER AND MOTHER ##

FATHER_BED_FINAL=$CONF_BED_DIR/${CROSS}_male_par.conf.1n_repeats_removed.bed
MOTHER_BED_FINAL=$CONF_BED_DIR/${CROSS}_fem_par.conf.1n_repeats_removed.bed

## FATHER
bedtools subtract -a $FATHER_BED -b $COMBINED_MASK > $FATHER_BED_FINAL

## MOTHER
bedtools subtract -a $MOTHER_BED -b $COMBINED_MASK > $MOTHER_BED_FINAL

#################################################################################
################ CREATE OFFSPRING CONFIDENT_REGIONS BED FILES ###################
#################################################################################

# Confident regions in the offspring should not be defined by the evidence in the offspring themsleves
# as we are not using any evidence in the offspring. Instead, these will also be defined by the parent's data

# We know where we are confident in the parents that there are no variants. Thus, excluding the extremely low likelihood of de novo variants, we can assume that these regions are also monomorphic in the offspring. 

# To find the offspring confident regions I will therefore take the intersection of the parents confident regions bed file. 

# E.g. the below confident region situation could occur

# MOTHER = Chr1:100-115
# FATHER = Chr1:107-122
#  Thus
# OFFSRPING = Chr1:107-115

# I will do this with bedtools intersect. 

## OFFSPRING
bedtools intersect -a $FATHER_BED -b $MOTHER_BED > $OFFSPRING_BED






#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --export=NONE
#SBATCH --array=1-5
#SBATCH --job-name=FINALISE_OFFSPRING_CONF_REGIONS
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

OFFSPRING_HOM_CONF_REGION_BED=$CONF_BED_DIR/${CROSS}_male.conf.bed

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

# >> In summary, to convert VCF coords -> BED coords, we -1 from the start point. << # 


#####################################################################################
######## Adding the variant positions to the offspring conf-regions bedfile #########
#####################################################################################

OFFSPRING_FILTERED_VARS=${WD}/Filtered_variants/${CROSS}_male_1.ALL_TRUTH_VARS.vcf.gz
OFFSPRING_VAR_POSITIONS_BED=${WD}/Filtered_variants/${CROSS}_var_positions.bed
OFFSPRING_CONF_WITH_VARS_BED=${CONF_BED_DIR}/${CROSS}_conf_regions_inc_vars.bed

## Get the positions of variants (1st two columns of vcf)
zcat $OFFSPRING_FILTERED_VARS | grep -v '#' | cut -f-2 > ${WD}/Filtered_variants/${CROSS}_var_positions.tmp

## copy the position column
cut -f2 ${WD}/Filtered_variants/${CROSS}_var_positions.tmp > ${WD}/Filtered_variants/${CROSS}_var_positions_col2.tmp

## combine the columns 
paste ${WD}/Filtered_variants/${CROSS}_var_positions.tmp ${WD}/Filtered_variants/${CROSS}_var_positions_col2.tmp >  ${WD}/Filtered_variants/${CROSS}_var_positions_prebed.tmp

## subtract 1 base from the start column
cat ${WD}/Filtered_variants/${CROSS}_var_positions_prebed.tmp | awk '{print $1 "\t" ($2 - 1) "\t" $3}' > $OFFSPRING_VAR_POSITIONS_BED

rm ${WD}/Filtered_variants/*${CROSS}*tmp

## combine the conf regions and var positions into single bed file.

cat $OFFSPRING_HOM_CONF_REGION_BED $OFFSPRING_VAR_POSITIONS_BED | bedtools sort | bedtools merge  > $OFFSPRING_CONF_WITH_VARS_BED




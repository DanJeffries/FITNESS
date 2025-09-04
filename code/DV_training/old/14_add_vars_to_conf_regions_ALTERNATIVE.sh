#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --export=NONE
#SBATCH --array=1
#SBATCH --job-name=FINALISE_OFFSPRING_CONF_REGIONS
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0
module load BEDTools/2.30.0-GCC-10.3.0

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training
RESEARCH_DIR=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training
CROSSES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/crossIDs.txt

FILTERED_GVCF_DIR=$RESEARCH_DIR/6_Confident_Var_sites/Filtered_VCFs
CONF_NONVAR_REGIONS=$RESEARCH_DIR/7_Confident_nonVar_sites
CONF_BED_DIR=$RESEARCH_DIR/8_Confident_regions

if [ ! -d "$CONF_BED_DIR" ]; then
   mkdir $CONF_BED_DIR
fi

CROSS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $CROSSES)

#####################################################################################
######## Adding the variant positions to the offspring conf-regions bedfile #########
#####################################################################################

CROSS_FILTERED_VARS=$FILTERED_GVCF_DIR/${CROSS}.ALL_TRUTH_VARS_CLEAN.vcf.gz
CROSS_VAR_POSITIONS_BED=${WD}/Filtered_variants/${CROSS}_var_positions.bed
CROSS_CONF_WITH_VARS_BED=${CONF_BED_DIR}/${CROSS}_conf_regions_inc_vars_ALT.bed

## Get the positions of variants (1st two columns of vcf)
zcat $CROSS_FILTERED_VARS | grep -v '#' | cut -f-2 > ${WD}/Filtered_variants/${CROSS}_var_positions.tmp

## copy the position column
cut -f2 ${WD}/Filtered_variants/${CROSS}_var_positions.tmp > ${WD}/Filtered_variants/${CROSS}_var_positions_col2.tmp

## combine the columns 
paste ${WD}/Filtered_variants/${CROSS}_var_positions.tmp ${WD}/Filtered_variants/${CROSS}_var_positions_col2.tmp >  ${WD}/Filtered_variants/${CROSS}_var_positions_prebed.tmp

## subtract 1 base from the start column
cat ${WD}/Filtered_variants/${CROSS}_var_positions_prebed.tmp | awk '{print $1 "\t" ($2 - 1) "\t" $3}' > $CROSS_VAR_POSITIONS_BED

#rm ${WD}/Filtered_variants/*${CROSS}*tmp

## combine the conf regions and var positions into single bed file.

CROSS_NONVAR_CONF_REGION_BED=$RESEARCH_DIR/7_Confident_nonVar_sites/Non_var_confident_regions/${CROSS}_CONF_HOM_REF_ALT_masked.bed


cat $CROSS_NONVAR_CONF_REGION_BED $CROSS_VAR_POSITIONS_BED | bedtools sort | bedtools merge  > $CROSS_CONF_WITH_VARS_BED

## Make a quick summary of the amount of sequence used for each sample

REF=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_noY.fna
SUMMARY=$RESEARCH_DIR/8_Confident_regions/${CROSS}_conf_regions_summary_ALT.txt

echo "$CROSS confident region summary" > $SUMMARY

N_conf_bases=$(bedtools getfasta -fi $REF -bed $CROSS_CONF_WITH_VARS_BED | grep -v '>' | wc -m)

echo "Total conf region bases (nonVar + Var): $N_conf_bases" >> $SUMMARY

N_non_var_bases=$(bedtools getfasta -fi $REF -bed $CROSS_NONVAR_CONF_REGION_BED | grep -v '>' | wc -m)

echo "Total conf nonVar region bases: $N_non_var_bases" >> $SUMMARY

N_var_bases=$(bedtools getfasta -fi $REF -bed $CROSS_VAR_POSITIONS_BED | grep -v '>' | wc -m)

echo "Total conf Var position bases: $N_var_bases" >> $SUMMARY




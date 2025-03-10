#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --array=1-5
#SBATCH --job-name=Conf_filter
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training
RESEARCH=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training
GVCFs=$RESEARCH/GVCF
MERGED_GVCFs=$WD/GVCFs_merged
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

#####################################################################################
## CONCAT chromosome GVCFs for each samaple
###########################################

FATHER_MERGED_GVCF=$MERGED_GVCFs/$FATHER.g.vcf.gz
MOTHER_MERGED_GVCF=$MERGED_GVCFs/$MOTHER.g.vcf.gz
OFFSPRING_MERGED_GVCF=$MERGED_GVCFs/$OFFSPRING.g.vcf.gz

bcftools concat -O z $GVCFs/$FATHER*chromosome*gz > $FATHER_MERGED_GVCF
bcftools concat -O z $GVCFs/$MOTHER*chromosome*gz > $MOTHER_MERGED_GVCF
bcftools concat -O z $GVCFs/$OFFSPRING*chromosome*gz > $OFFSPRING_MERGED_GVCF

FATHERS_CONF_GVCF=$FILTERED_GVCF_DIR/$FATHER.conf.g.vcf.gz
MOTHERS_CONF_GVCF=$FILTERED_GVCF_DIR/$MOTHER.conf.g.vcf.gz
OFFSPRING_CONF_GVCF=$FILTERED_GVCF_DIR/$OFFSPRING.conf.g.vcf.gz

#####################################################################################
## Filters to keep confident regions 
#####################################

# homozygous
# minGQ (correct genotype) = 30
# min depth = 15

## Filter for father
bcftools view $FATHER_MERGED_GVCF -i 'GT="hom"' | \
bcftools view -i 'FORMAT/DP >= 15' $FATHER_GVCF | \
bcftools view -i 'FORMAT/GQ >= 30' -O z > $FATHERS_CONF_GVCF

## Filter for Mother
bcftools view $MOTHER_MERGED_GVCF -i 'GT="hom"' | \
bcftools view -i 'FORMAT/DP >= 15' $MOTHER_GVCF | \
bcftools view -i 'FORMAT/GQ >= 30' -O z > $MOTHERS_CONF_GVCF


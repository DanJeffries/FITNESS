#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --array=1-5
#SBATCH --job-name=BCFtools_VarFilter
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training
STORAGE=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training
CROSSES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/crossIDs.txt
JOINT_CALL_VCF=$STORAGE/4_GATK/Unfiltered_VCF/Joint_all_unfiltered.vcf.gz
FILTERED_VCF_DIR=$WD/Filtered_variants

if [ ! -d "$FILTERED_VCF_DIR" ]; then
   mkdir $FILTERED_VCF_DIR
fi

CROSS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $CROSSES)

FATHER="${CROSS}_male_par"
MOTHER="${CROSS}_fem_par"

FATHERS_VCF=$FILTERED_VCF_DIR/$FATHER.vcf.gz
MOTHERS_VCF=$FILTERED_VCF_DIR/$MOTHER.vcf.gz


## Filters (required to KEEP loci)

## GATK hard filters (https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants)
# FS (Fishers exact test for strand bias) < 60
# QD (Variant confidence divided by unfiltered depth) > 2.0 
# SOR (Strand Odds Ratio) < 3.0 
# MQ (Root mean square of mapping quality for all reads at the site) > 40
# MQRankSum (Compares mapping quality of reads supporting the REF vs ALT allele) > -12.5
# ReadPosRankSum (position of variants in reads) > -8.0

## In addition to the GATK recomended hard filters above we also need to filter by the following genotype confidence related stats. 

# Both parents homozygous
# Biallelic 
# minQ (polymorphism exists) = 30
# minGQ (correct genotype) = 30
# min depth = 20
# zero reads in support of non-called allele. 
#

## Note that we cannot subset by sample and filter at the same time, because the filters are done before the sample subsetting in the order of operations. But we can just pipe then! To make this more straight forward we will do this separately
## for each sample, making a separate VCF for each.

## Filter for father
bcftools view -M 2 -s $FATHER $JOINT_CALL_VCF | \
bcftools filter -i 'INFO/QD >= 2.0 & INFO/FS <= 60 & INFO/SOR <= 3.0 & INFO/MQ >= 40.0 & INFO/MQRankSum >= -12.5 & INFO/ReadPosRankSum >= -8.0' | \
bcftools filter -i 'GT="hom"' | \
bcftools filter -i 'QUAL >= 30' | \
bcftools filter -i 'FORMAT/DP >= 20' | \
bcftools filter -i 'FORMAT/DP <= 100' | \
bcftools filter -i 'FORMAT/GQ >= 30' | \
bcftools filter -i 'FORMAT/AD[0:0]=0 | FORMAT/AD[0:1]=0' \
                -O z >  $FATHERS_VCF

tabix $FATHERS_VCF

## Filter for Mother
bcftools view -M 2 -s $MOTHER $JOINT_CALL_VCF | \
bcftools filter -i 'INFO/QD >= 2.0 & INFO/FS <= 60 & INFO/SOR <= 3.0 & INFO/MQ >= 40.0 & INFO/MQRankSum >= -12.5 & INFO/ReadPosRankSum >= -8.0' | \
bcftools filter -i 'GT="hom"' | \
bcftools filter -i 'QUAL >= 30' | \
bcftools filter -i 'FORMAT/DP >= 20' | \
bcftools filter -i 'FORMAT/DP <= 100' | \
bcftools filter -i 'FORMAT/GQ >= 30' | \
bcftools filter -i 'FORMAT/AD[0:0]=0 | FORMAT/AD[0:1]=0' \
                -O z >  $MOTHERS_VCF

tabix $MOTHERS_VCF

# old filter - replaced by zero reads for non-called allele. 
#bcftools filter -i 'FORMAT/AD[0:0]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.85 | FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.85' \


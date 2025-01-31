#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --job-name=BCFtools_GATK_baseline
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training
STORAGE=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training
JOINT_CALL_VCF=$STORAGE/4_GATK/Unfiltered_VCF/Joint_all_unfiltered.vcf.gz
FILTERED_VCF_DIR=$STORAGE/10_Validation/GATK_Filtered_VCFs/gatk_baseline/

if [ ! -d "$FILTERED_VCF_DIR" ]; then
   mkdir $FILTERED_VCF_DIR
fi

OFFSPRING="FG_male_1"
OFFSPRING_VCF=$FILTERED_VCF_DIR/${OFFSPRING}.GATK.filtered_baseline.vcf.gz

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

TEST_REGIONS=$WD/training_regions/FG_test_partitions.bed

## Filter for father
bcftools view   -M 2 \
                -s $OFFSPRING       \
		-R $TEST_REGIONS    \
		   $JOINT_CALL_VCF| \
bcftools filter -i 'INFO/QD >= 2.0 & INFO/FS <= 60 & INFO/SOR <= 3.0 & INFO/MQ >= 40.0 & INFO/MQRankSum >= -12.5 & INFO/ReadPosRankSum >= -8.0' | \
bcftools filter -i 'QUAL >= 30'       | \
bcftools filter -i 'FORMAT/DP >= 20'  | \
bcftools filter -i 'FORMAT/DP <= 100' | \
bcftools filter -i 'FORMAT/GQ >= 30'    \
                -O z >  $OFFSPRING_VCF

tabix $OFFSPRING_VCF

#bcftools filter -i 'FORMAT/AD[0:0]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.85 | FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.85' \


#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --job-name=BCFtools_Filter
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load vital-it
module load UHTS/Analysis/vcftools/0.1.15

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/Find_2n_X_windows/
JOINT_CALL_VCF=$WD/Unfiltered_VCF/Joint_all_unfiltered.vcf
FILTERED_VCF_DIR=$WD/Filtered_VCF

if [ ! -d "$FILTERED_VCF_DIR" ]; then
   mkdir $FILTERED_VCF_DIR
fi


## Filters (required to KEEP loci)

## GATK hard filters (https://gatk.broadinstitute.org/hc/en-us/articles/360035890471-Hard-filtering-germline-short-variants)
# FS (Fishers exact test for strand bias) < 60
# QD (Variant confidence divided by unfiltered depth) > 2.0 
# SOR (Strand Odds Ratio) < 3.0 
# MQ (Root mean square of mapping quality for all reads at the site) > 40
# MQRankSum (Compares mapping quality of reads supporting the REF vs ALT allele) > -12.5
# ReadPosRankSum (position of variants in reads) > -8.0

## In addition to the GATK recomended hard filters above we also need to filter by the following genotype confidence related stats. 

# -M (max alleles = 2) = 2
# minQ (polymorphism exists) = 30
# minGQ (correct genotype) = 30
# min depth = 10
# 0.15 <= allele balance <= 0.85

## Note that we cannot subset by sample and filter at the same time, because the filters are done before the sample subsetting in the order of operations. But we can just pipe then! To make this more straight forward we will do this separately
## for each sample, making a separate VCF for each.

## Filter 
bcftools view -M 2 $JOINT_CALL_VCF | \
bcftools filter -i 'INFO/QD >= 2.0 & INFO/FS <= 60 & INFO/SOR <= 3.0 & INFO/MQ >= 40.0 & INFO/MQRankSum >= -12.5 & INFO/ReadPosRankSum >= -8.0' | \
bcftools filter -i 'QUAL >= 30' | \
bcftools filter -i 'FORMAT/AD[0:0] + FORMAT/AD[0:1] >= 10' | \
bcftools filter -i 'FORMAT/GQ >= 30' | \
bcftools filter -i 'FORMAT/AD[0:0]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.85 | FORMAT/AD[0:1]/(FORMAT/AD[0:0]+FORMAT/AD[0:1]) >= 0.85' >  $FILTERED_VCF_DIR/Filtered_BiAL_QUAL_GQ30_MinDP10_AF085_HARD.vcf




#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=02:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem=110G
#SBATCH --export=NONE
#SBATCH --array=1-5
#SBATCH --job-name=IMPUTE_OFFSPRING_GTs
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/
CROSSES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/crossIDs.txt
CROSS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $CROSSES)

FATHERS_FILTERED_VARIANTS=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/Filtered_variants/${CROSS}_male_par.vcf.gz
MOTHERS_FILTERED_VARIANTS=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/Filtered_variants/${CROSS}_fem_par.vcf.gz

#########
## Step 1 - Get HOM_REF and HOM_ALT locus lists for each parent. 
##################################################################################################

FATHER_HOM_REF=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/Filtered_variants/${CROSS}_male_par_HOM_REF
FATHER_HOM_ALT=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/Filtered_variants/${CROSS}_male_par_HOM_ALT

#vcftools --gzvcf $FATHERS_FILTERED_VARIANTS --max-non-ref-ac 0 --counts2 --out $FATHER_HOM_REF
#vcftools --gzvcf $FATHERS_FILTERED_VARIANTS --non-ref-ac 2 --counts2 --out $FATHER_HOM_ALT

MOTHER_HOM_REF=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/Filtered_variants/${CROSS}_fem_par_HOM_REF
MOTHER_HOM_ALT=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/Filtered_variants/${CROSS}_fem_par_HOM_ALT

#vcftools --gzvcf $MOTHERS_FILTERED_VARIANTS --max-non-ref-ac 0 --counts2 --out $MOTHER_HOM_REF
#vcftools --gzvcf $MOTHERS_FILTERED_VARIANTS --non-ref-ac 2 --counts2 --out $MOTHER_HOM_ALT

########
## Step 2 - Find the loci with HOM_REF in one parent and HOM_ALT in another, or HOM_ALT in both
##################################################################################################

OFFSRPING_HET_LIST=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/Filtered_variants/${CROSS}_male_HET.list
OFFSRPING_HOM_ALT_LIST=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/Filtered_variants/${CROSS}_male_HOM_ALT.list

## FATHER: HOM_REF, MOTHER: HOM_ALT
#comm -12 <(cut -f1,2 $FATHER_HOM_REF.frq.count | grep -v 'CHROM'|sort) <(cut -f1,2 $MOTHER_HOM_ALT.frq.count | grep -v 'CHROM'| sort) > $OFFSRPING_HET_LIST

## FATHER: HOM_ALT, MOTHER: HOM_REF
#comm -12 <(cut -f1,2 $FATHER_HOM_ALT.frq.count | grep -v 'CHROM'|sort) <(cut -f1,2 $MOTHER_HOM_REF.frq.count | grep -v 'CHROM'| sort) >> $OFFSRPING_HET_LIST

## FATHER: HOM_ALT, MOTHER: HOM_ALT
#comm -12 <(cut -f1,2 $FATHER_HOM_ALT.frq.count | grep -v 'CHROM'|sort) <(cut -f1,2 $MOTHER_HOM_ALT.frq.count | grep -v 'CHROM'| sort) > $OFFSRPING_HOM_ALT_LIST

#####
## Step 3 - Filter the offspring GVCF for the kept loci 
##################################################################################################

# Note we do this separately for the HET and HOM_ALT, this way I can more easily change the genotypes in the offspring's VCF

GVCFs=$WD/GVCF
OFFSPRING_MERGED_GVCF=$WD/GVCFs_merged/${CROSS}_male_1.g.vcf.gz

#bcftools concat -O z $GVCFs/${CROSS}_male_1_*chromosome*gz > $OFFSPRING_MERGED_GVCF

OFFSPRING_HET_VCF=$WD/Filtered_variants/${CROSS}_male_1.HET.vcf
OFFSPRING_HOM_ALT_VCF=$WD/Filtered_variants/${CROSS}_male_1.HOM_ALT.vcf

#vcftools --gzvcf $OFFSPRING_MERGED_GVCF --positions $OFFSRPING_HET_LIST --recode --out $OFFSPRING_HET_VCF
#vcftools --gzvcf $OFFSPRING_MERGED_GVCF --positions $OFFSRPING_HOM_ALT_LIST --recode --out $OFFSPRING_HOM_ALT_VCF

######
## Step 4. Change any loci that don't match the expectation based on the parents. 
#########################################################################################################################

## HET
# Change any 0/0, 0|0, 1/1, 1|1 to 0/1

OFFSPRING_HET_CORRECTED_VCF=$WD/Filtered_variants/${CROSS}_male_1.HET.CORRECTED.vcf

#sed 's/0\/0/0\/1/g' $OFFSPRING_HET_VCF.recode.vcf | sed 's/0|0/0|1/g' | sed 's/1\/1/0\/1/g' | sed 's/1|1/0|1/g' > $OFFSPRING_HET_CORRECTED_VCF

## HOM_ALT
# Change any 0/0, 0|0, 0/1, 0|1 to 1/1

OFFSPRING_HOM_ALT_CORRECTED_VCF=$WD/Filtered_variants/${CROSS}_male_1.HOM_ALT.CORRECTED.vcf

#sed 's/0\/0/1\/1/g' $OFFSPRING_HOM_ALT_VCF.recode.vcf | sed 's/0|0/1|1/g' | sed 's/0\/1/1\/1/g' | sed 's/0|1/1|1/g' > $OFFSPRING_HOM_ALT_CORRECTED_VCF


## Step 5. Concatenate the HET and HOM_ALT corrected VCFs, gzip and index them. 

OFFSPRING_ALL_TRUTH_VARS_VCF=$WD/Filtered_variants/${CROSS}_male_1.ALL_TRUTH_VARS.CORRECTED.vcf

bcftools concat $OFFSPRING_HET_CORRECTED_VCF $OFFSPRING_HOM_ALT_CORRECTED_VCF > $OFFSPRING_ALL_TRUTH_VARS_VCF

bgzip $OFFSPRING_ALL_TRUTH_VARS_VCF

tabix $OFFSPRING_ALL_TRUTH_VARS_VCF.gz





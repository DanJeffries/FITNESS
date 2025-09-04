#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=01:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=50G
#SBATCH --export=NONE
#SBATCH --array=1-5
#SBATCH --job-name=IDENTIFY_TRUTH_VARIANTS
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

## Here we want to identify true variants. We will do this by identifying loci which are fixed for different alleles in the parents (e.g. mother = 0/0, father = 1/1) and where the offspring's genotype matches the expectation (e.g. 0/1). 

module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/
CROSSES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/crossIDs.txt
CROSS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $CROSSES)

FATHERS_FILTERED_VARIANTS=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/Filtered_variants/${CROSS}_male_par.vcf
MOTHERS_FILTERED_VARIANTS=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/Filtered_variants/${CROSS}_fem_par.vcf

#########
## Step 1 - Get HOM_REF and HOM_ALT locus lists for each parent. 
##################################################################################################

FATHER_HOM_REF=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/Filtered_variants/${CROSS}_male_par_HOM_REF
FATHER_HOM_ALT=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/Filtered_variants/${CROSS}_male_par_HOM_ALT

#vcftools --gzvcf $FATHERS_FILTERED_VARIANTS.gz --max-non-ref-ac 0 --recode --out $FATHER_HOM_REF
#vcftools --gzvcf $FATHERS_FILTERED_VARIANTS.gz --non-ref-ac 2 --recode --out $FATHER_HOM_ALT

MOTHER_HOM_REF=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/Filtered_variants/${CROSS}_fem_par_HOM_REF
MOTHER_HOM_ALT=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/Filtered_variants/${CROSS}_fem_par_HOM_ALT

#vcftools --gzvcf $MOTHERS_FILTERED_VARIANTS.gz --max-non-ref-ac 0 --recode --out $MOTHER_HOM_REF
#vcftools --gzvcf $MOTHERS_FILTERED_VARIANTS.gz --non-ref-ac 2 --recode --out $MOTHER_HOM_ALT

#######s 
## Step 2 - Find the loci with HOM_REF in one parent and HOM_ALT in another, or HOM_ALT in both
##################################################################################################

OFFSRPING_HET_LIST=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/Filtered_variants/${CROSS}_male_1_HET.list
OFFSRPING_HOM_ALT_LIST=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/Filtered_variants/${CROSS}_male_1_HOM_ALT.list

## FATHER: HOM_REF, MOTHER: HOM_ALT
#comm -12 <(grep -v '#' ${FATHER_HOM_REF}.recode.vcf |  cut -f1,2 | sort ) <(grep -v '#' ${MOTHER_HOM_ALT}.recode.vcf |  cut -f1,2 | sort) > $OFFSRPING_HET_LIST

## FATHER: HOM_ALT, MOTHER: HOM_REF
#comm -12 <(grep -v '#' ${FATHER_HOM_ALT}.recode.vcf | cut -f1,2 | sort) <(grep -v '#' ${MOTHER_HOM_REF}.recode.vcf | cut -f1,2 | sort) >> $OFFSRPING_HET_LIST

## FATHER: HOM_ALT, MOTHER: HOM_ALT
#comm -12 <(grep -v '#' ${FATHER_HOM_ALT}.recode.vcf | cut -f1,2 | sort) <(grep -v '#' ${MOTHER_HOM_ALT}.recode.vcf | cut -f1,2 | sort) > $OFFSRPING_HOM_ALT_LIST

## bgzip and index the files now we are done with them

#bgzip $FATHER_HOM_REF.recode.vcf
#bgzip $FATHER_HOM_ALT.recode.vcf

#tabix $FATHER_HOM_REF.recode.vcf.gz
#tabix $FATHER_HOM_ALT.recode.vcf.gz

#bgzip $MOTHER_HOM_REF.recode.vcf
#bgzip $MOTHER_HOM_ALT.recode.vcf

#tabix $MOTHER_HOM_REF.recode.vcf.gz
#tabix $MOTHER_HOM_ALT.recode.vcf.gz


#####
## Step 3 - Create individual sample VCFs from the Joint called VCF with all samples. For each sample, keep only loci identified for that cross in the list above as either fixed homozygous in parents or both parents being HOM_ALT. These will either be HET or HOM_ALT in the offspring. 
##################################################################################################

# Note we do this separately for the HET and HOM_ALT, this way I can more easily filter for only the genotypes in the offspring's VCF that match mendelian expectations

JOINT_VCF=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training/4_GATK/Unfiltered_VCF/Joint_all_unfiltered.vcf.gz

OFFSPRING_HET_VCF=$WD/Filtered_variants/${CROSS}_male_1.HET.vcf.gz
OFFSPRING_HOM_ALT_VCF=$WD/Filtered_variants/${CROSS}_male_1.HOM_ALT.vcf.gz

### first for het positions in offspring
#bcftools view -s ${CROSS}_male_1 \
#	      -R $OFFSRPING_HET_LIST \
#	         $JOINT_VCF |\
#bcftools view -g het \
#              -O z \
#              > $OFFSPRING_HET_VCF ## pipe to this second command to keep only loci that agree with mendelian expectations. Don't do this all in first command because GT filters MAY get done before the sample filtering. 

#tabix $OFFSPRING_HET_VCF

## now for hom_alt
#bcftools view -s ${CROSS}_male_1 \
#              -R $OFFSRPING_HOM_ALT_LIST \
#                 $JOINT_VCF |\
#bcftools view -g hom \
#              -O z \
#              > $OFFSPRING_HOM_ALT_VCF
#tabix $OFFSPRING_HOM_ALT_VCF

## Step 5. Concatenate the HET and HOM_ALT corrected VCFs, gzip and index them. 

OFFSPRING_ALL_TRUTH_VARS_VCF=$WD/Filtered_variants/${CROSS}_male_1.ALL_TRUTH_VARS.vcf.gz

#bcftools concat -a $OFFSPRING_HET_VCF $OFFSPRING_HOM_ALT_VCF \
#	        -O u | \
#bcftools view   -a \
#		-O z \
#		> $OFFSPRING_ALL_TRUTH_VARS_VCF

#tabix $OFFSPRING_ALL_TRUTH_VARS_VCF

## Step 6. Make a VCF for the cross for validation / summary plots. 

ALL_TRUTH_VARS_LIST=$WD/Filtered_variants/${CROSS}.ALL_TRUTH_VARS.list

#zcat $OFFSPRING_ALL_TRUTH_VARS_VCF | grep -v '#' | cut -f1,2 > $ALL_TRUTH_VARS_LIST

CROSS_TRUTH_VAR_VCF=$WD/Filtered_variants/${CROSS}.ALL_TRUTH_VARS_CLEAN.vcf.gz

echo "${CROSS}_male_par,${CROSS}_fem_par,${CROSS}_male_1" 

#bcftools view $JOINT_VCF \
#              -s ${CROSS}_male_par,${CROSS}_fem_par,${CROSS}_male_1 \
#              -R $ALL_TRUTH_VARS_LIST \
#	      -O z \
#              > $CROSS_TRUTH_VAR_VCF \

#tabix $CROSS_TRUTH_VAR_VCF

## Step 7. Make a VCF just for the offspring, which will be the input for the training

OFFSPRING_TRUTH_VAR_VCF=$WD/Filtered_variants/${CROSS}_male_1.ALL_TRUTH_VARS_CLEAN.SNPS_ONLY.vcf.gz

bcftools view $CROSS_TRUTH_VAR_VCF \
	      -s ${CROSS}_male_1 \
	      -v snps \
	      -O z \
	      > $OFFSPRING_TRUTH_VAR_VCF


#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --array=1-5
#SBATCH --job-name=Make_non_var_conf_regions_pedigree_GVCF
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0

bcftools=/storage/homefs/dj20y461/Software/bcftools-1.21/bcftools

## The purpose of this script is to identify non-variant confident regions in the offspring based on the evidence in the parents. The logic is as follows:

# In order to train the model to learn what a non-variant position looks like, the "confident regions" must contain loci which we are confident are non-variant in the offspring. However, we do not want to use
# evidence in the offspring to identify these positions. For example, if we removed loci called as 0/1 in the offspring (despite parents both being 0/0) we would be removing false positive variant signal. We want DV to learn what false positives
# look like, so we want these situations in there. 

# The solution is to identify genomic regions that we EXPECT to be homozygous in the offspring based on the parents genotypes. Such loci should have identical genotypes in the parents, which are called with extremely good confidence. 
# Note that it is not enough just that both parents be homozygous, they must be homozygous for the same allele. E.g. AA & AA, not AA & TT. 

# The below filters will identify such regions in the parents:

###############################
## >>> Set paths & vars  <<< ##
###############################

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training
MERGED_GVCFs=$WD/GVCFs_merged
CONFIDENT_REGIONS=$WD/Confident_regions
CROSSES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/crossIDs.txt

CROSS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $CROSSES)

###################################################
### >>> Combine parent and offspring filtered GVCFs
##################################################

PARENTS_GVCF=$MERGED_GVCFs/${CROSS}_parents.g.vcf.gz
CROSS_CONF_REGIONS_GVCF=$CONFIDENT_REGIONS/${CROSS}_male_1_putative_nonVar_confident_regions.g.vcf.gz
CROSS_CONF_BED_FINAL=$CONFIDENT_REGIONS/${CROSS}_nonVar_conf_regions_masked.bed

PEDIGREE_nonVAR=$CONFIDENT_REGIONS/${CROSS}_putative_nonVar.g.vcf.gz
PEDIGREE_CONF_REGIONS_GVCF=$CONFIDENT_REGIONS/${CROSS}_putative_nonVar_confident_regions.g.vcf.gz

#bcftools merge $PARENTS_GVCF \
#	       $CROSS_CONF_REGIONS_GVCF \
#	       -O z \
#	       > $PEDIGREE_nonVAR 

#tabix $PEDIGREE_nonVAR

bcftools view  $PEDIGREE_nonVAR \
               -R $CROSS_CONF_BED_FINAL \
	       -i "N_MISSING=0" \
               -a \
               -O z \
               > $PEDIGREE_CONF_REGIONS_GVCF

#tabix $PEDIGREE_CONF_REGIONS_GVCF

### do i need to add a missing values filter in? I think not, DV takes the bam file, so if there is no call in the offspring, it means that there were not enough reads for GATK, but there might be for DV. Eitherway, it is a real world representation of what homozygous regions can look like.





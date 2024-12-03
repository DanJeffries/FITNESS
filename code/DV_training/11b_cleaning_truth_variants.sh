#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=01:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=5
#SBATCH --mem=50G
#SBATCH --export=NONE
#SBATCH --array=1-5
#SBATCH --job-name=CLAEN_TRUTH_VARIANTS
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

## Despite the filtering already done there are still a handful of loci (a few hundred per cross) that aren't suitable for training. They probably wouldn't make a big difference, but for the sake of completeness I will remove them anyway. 

module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training/
CROSSES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/crossIDs.txt
CROSS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $CROSSES)


CROSS_TRUTH_VAR_VCF=$WD/Filtered_variants/${CROSS}.ALL_TRUTH_VARS.vcf.gz
FAILS=$WD/Filtered_variants/${CROSS}_fails.tsv

CROSS_TRUTH_VAR_VCF_CLEAN=$WD/Filtered_variants/${CROSS}.ALL_TRUTH_VARS_CLEAN.vcf.gz

bcftools view $CROSS_TRUTH_VAR_VCF \
              -T ^${FAILS} \
	      -O z \
              > $CROSS_TRUTH_VAR_VCF_CLEAN \

tabix $CROSS_TRUTH_VAR_VCF_CLEAN



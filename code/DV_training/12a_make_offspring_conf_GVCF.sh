#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=06:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --export=NONE
#SBATCH --array=1-5
#SBATCH --job-name=Make_offspring_conf_regions_GVCF
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0
module load BEDTools/2.30.0-GCC-10.3.0

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training
CROSSES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/crossIDs.txt

CONF_REGIONS_DIR=$WD/Confident_regions

CROSS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $CROSSES)
SAMPLE=${CROSS}_male_1

MERGED_GVCF=$WD/GVCFs_merged/${SAMPLE}.g.vcf.gz
OFFSPRING_NONVAR_CONF_GVCF=$CONF_REGIONS_DIR/${SAMPLE}_putative_nonVar_confident_regions.g.vcf.gz

### I want to see what the confident regions resulting from my filters look like. Specifically how many things that could look like SNPs are in there? 
### What is the depth distribution like
### What are the GQ distributions like?
### How many things that look like SNPs to GATK are there in there? 
### What are the ADR values like

## I have made the offspring non_var confident region gVCF, so now I would like to generate some stats for this

bcftools stats $OFFSPRING_NONVAR_CONF_GVCF  > $CONF_REGIONS_DIR/${SAMPLE}_putative_nonVar_confident_regions.stats






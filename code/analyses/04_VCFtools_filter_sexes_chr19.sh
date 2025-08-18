#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=04:00:00 
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --export=NONE
##SBATCH --array=1-22
#SBATCH --job-name=VCFtools_FILTER_females_chr19
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load VCFtools/0.1.16-GCC-10.3.0

WD=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS

SEX="female"
VCF=$WD/DV_calling/CHR_VCFs_prefiltered_NOFILL/chr_19.allsamples.prefiltered.vcf.gz
POPMAP=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/sample_metadata/inferred_sexes/females.txt

## make the output directory
FILTER_OUTDIR=$WD/DV_calling/chr19_sepsexes/$SEX

if [ ! -d "$FILTER_OUTDIR" ]; then
	mkdir -p $FILTER_OUTDIR
fi

## calculate statistics


vcftools --gzvcf $VCF \
	 --keep $POPMAP \
	 --maf 0.05 \
	 --remove-indels \
	 --max-missing 0.2 \
	 --recode \
	 --out $FILTER_OUTDIR/all_samples_chr_19_${SEX}s


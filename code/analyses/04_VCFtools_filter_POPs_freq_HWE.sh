#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=04:00:00 
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=64G
#SBATCH --export=NONE
#SBATCH --array=1-22
#SBATCH --job-name=VCFtools_FILTER_FREQ_HWE
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load VCFtools/0.1.16-GCC-10.3.0

WD=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS

#FG  LG  SL  SR  TL  WB  WK  WT
POP=$1

VCF=$WD/DV_calling/CHR_VCFs_prefiltered_NOFILL/chr_${SLURM_ARRAY_TASK_ID}.allsamples.prefiltered.vcf.gz
POPMAP=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/analyses/popmaps/full/${POP}_samples.txt


## potential stats to use are below, can swap them out of the command as needed (only one per command is allowed)

#        --het \
#        --hardy \
#        --freq2
#        --SNPdensity \
#	 --missing-indv \
#        --missing-site \
#        --site-mean-depth \
#	 --extract-FORMAT-info AD \

## make the output directory
FILTER_OUTDIR=$WD/DV_calling/POP_FILTERED_FREQ_HWE/$POP

if [ ! -d "$FILTER_OUTDIR" ]; then
	mkdir -p $FILTER_OUTDIR
fi

## calculate statistics


vcftools --gzvcf $VCF \
	 --keep $POPMAP \
	 --maf 0.1 \
	 --remove-indels \
	 --hwe 0.05 \
	 --max-missing 0.2 \
	 --recode \
	 --out $FILTER_OUTDIR/${POP}_chr_${SLURM_ARRAY_TASK_ID}.filtered_freq_hwe


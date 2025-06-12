#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=01:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=4G
#SBATCH --export=NONE
#SBATCH --job-name=CONCAT_VCFs_cohort_1
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

#### script to concatenate VCFs by interval

module load BCFtools/1.12-GCC-10.3.0

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_calling
COHORT=cohort_1
VCFs=$WD/GLnexus/$COHORT/*vcf.gz

OUTDIR=$WD/concatenated_cohort_VCFs

if [ ! -d "$OUTDIR" ]; then
        mkdir -p $OUTDIR
fi

OUT_VCF=$OUTDIR/${COHORT}_concatenated.vcf.gz

N_THREADS=40

bcftools concat \
	 $VCFs \
	 -O z \
	 --threads $N_THREADS \
	 > $OUT_VCF

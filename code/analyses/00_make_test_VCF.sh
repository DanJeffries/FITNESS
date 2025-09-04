#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=23:59:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --job-name=subsample
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load Anaconda3
eval "$(conda shell.bash hook)"

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS

GZVCF=$WD/DV_calling/concatenated_cohort_VCFs/cohort_1_concatenated.vcf.gz

r=0.01

SUBSAMPLED_GZVCF=$WD/analyses/cohort_1_tests/cohort_1_concatenated.subsampled_${r}.vcf.gz

conda activate vcflib_env

bcftools view $GZVCF |\
\
vcfrandomsample -r $r |\
\
bcftools view -O z > $SUBSAMPLED_GZVCF 



#usage: vcfrandomsample [options] [<vcf file>]

#options:
#    -r, --rate RATE          base sampling probability per locus
#    -s, --scale-by KEY       scale sampling likelihood by this Float info field
#    -p, --random-seed N      use this random seed (by default read from /dev/random)
#    -q, --pseudorandom-seed  use a pseudorandom seed (by default read from /dev/random)
#
#Randomly sample sites from an input VCF file, which may be provided as stdin.
#Scale the sampling probability by the field specified in KEY.  This may be
#used to provide uniform sampling across allele frequencies, for instance.




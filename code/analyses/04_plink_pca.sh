#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=02:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem=99G
#SBATCH --export=NONE
#SBATCH --array=1
#SBATCH --job-name=plink_PCA
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

PLINK=~/Software/plink/plink

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS

MIXED_GZVCF=$WD/analyses/cohort_1_tests/mixed_cohort_1.subsample_0.01_filtered.recode.vcf.gz
OUT_PREFIX=$WD/analyses/cohort_1_tests/stats/cohort_1_downsampled

$PLINK --vcf $MIXED_GZVCF \
       --double-id \
       --allow-extra-chr \
       --set-missing-var-ids @:# \
       --make-bed \
       --pca \
       --vcf-half-call m \
       --mind \
       --out $OUT_PREFIX


#  plink <input flag(s)...> [command flag(s)...] [other flag(s)...]
#  plink --help [flag name(s)...]
#
#Commands include --make-bed, --recode, --flip-scan, --merge-list,
#--write-snplist, --list-duplicate-vars, --freqx, --missing, --test-mishap,
#--hardy, --mendel, --ibc, --impute-sex, --indep-pairphase, --r2, --show-tags,
#--blocks, --distance, --genome, --homozyg, --make-rel, --make-grm-gz,
#--rel-cutoff, --cluster, --pca, --neighbour, --ibs-test, --regress-distance,
#--model, --bd, --gxe, --logistic, --dosage, --lasso, --test-missing,
#--make-perm-pheno, --tdt, --qfam, --annotate, --clump, --gene-report,
#--meta-analysis, --epistasis, --fast-epistasis, and --score.

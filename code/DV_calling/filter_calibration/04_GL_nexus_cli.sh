#!/bin/bash

#SBATCH --partition=epyc2
#SBATCH --time=8:00:00
#SBATCH --tasks=1
#SBATCH --cpus-per-task=20
#SBATCH --mem-per-cpu=67G
##SBATCH --mem=100G
#SBATCH --export=NONE
#SBATCH --array=1
#SBATCH --job-name=GLNEXUS_cohort_FG_family
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load BCFtools/1.12-GCC-10.3.0
module load HTSlib/1.12-GCC-10.3.0

#### Joint calling with GLnexus
GLNEXUS=~/Software/glnexus_cli  ## using precompiled executable v1.4.1 from here: https://github.com/dnanexus-rnd/GLnexus/releases 

#WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_calling

WD=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_calling/filter_calibration/

################################
#####>>>>>> COHORT  <<<<<#######
################################

COHORT=FG

## Set directories

POSTPROCESS_DIR=$WD/postprocess_variants/${COHORT}* ## input dir 
JOINT_CALLING_OUTDIR=$WD/GLnexus ## output dir

## splitting the joint calling into cohorts and by chromosome
CHROMOSOME_BED=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_calling/chr_beds/chr_${SLURM_ARRAY_TASK_ID}.bed
CHROMOSOME=$(cut -f1 $CHROMOSOME_BED)

## Fix the last chrom name, as this is a bunch of unplaced contigs. 
if (($SLURM_ARRAY_TASK_ID == 22)); then
	CHROMOSOME="unplaced_contigs"
fi 

## make separate output dir for each cohort
COHORT_OUTDIR=$JOINT_CALLING_OUTDIR/$COHORT ## output dir

if [ ! -d "$COHORT_OUTDIR" ]; then
        mkdir -p $COHORT_OUTDIR
fi

## Glnexus will fail if the temporary DB directory that it creates already exists. So make a unique one for each run - can remove after

GLnexus_DB_temp_dir=$JOINT_CALLING_OUTDIR/GL_nexus_tempDBs
GL_temp_DB=$GLnexus_DB_temp_dir/GLnexus.$COHORT.chr_${CHROMOSOME}.DB


if [ ! -d "$GLnexus_DB_temp_dir" ]; then
        mkdir -p $GLnexus_DB_temp_dir
fi

if [ -d "$GL_temp_DB" ]; then
        rm $GL_temp_DB* -rf
fi

## Set outfile path
OUT_PREFIX=$JOINT_CALLING_OUTDIR/${COHORT}/${COHORT}_chr_${CHROMOSOME}

echo "$GLNEXUS --config DeepVariant \
            $POSTPROCESS_DIR/*g.vcf.gz \
            --dir $GL_temp_DB \
            --bed $CHROMOSOME_BED \
            --mem-gbytes 1000 \
	    > $OUT_PREFIX.bcf"

## Run
$GLNEXUS --config DeepVariant \
            $POSTPROCESS_DIR/*g.vcf.gz \
	    --dir $GL_temp_DB \
	    --bed $CHROMOSOME_BED \
	    --mem-gbytes 1000 \
    	    > $OUT_PREFIX.bcf

## convert to gzvcf and index. 
bcftools view $OUT_PREFIX.bcf | bgzip -@ 20 -c > $OUT_PREFIX.vcf.gz
tabix $OUT_PREFIX.vcf.gz

## clean up temp dirs if needed. 
rm $GL_temp_DB* -rf



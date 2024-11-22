#!/bin/bash

#SBATCH --partition=bdw
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --tasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=32G
#SBATCH --export=NONE
#SBATCH --job-name=Mendelian_eval
#SBATCH --array=1-5
#SBATCH --output=%x_%A-%a.out
#SBATCH --error=%x_%A-%a.err

module load BCFtools/1.12-GCC-10.3.0
module load VCFtools/0.1.16-GCC-10.3.0
module load ant/1.10.11-Java-11 ## needed for RTG

## set path for RealTimeGenomics toolkit (installed from github: https://github.com/RealTimeGenomics/rtg-tools) 

RTG=~/Software/rtg-tools/dist/rtg-tools-3.12.1-32f94bb2/rtg

WD=/storage/scratch/iee/dj20y461/Stickleback/G_aculeatus/FITNESS/DV_training
CROSSES=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/crossIDs.txt

CROSS=$(sed -n "${SLURM_ARRAY_TASK_ID}p" < $CROSSES)

###################################################################
### Prepare ref fasta for RTG by converting to SDF
##################################################

REF=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames.fna
SDF_OUT=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/ref/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_SDF

if [ ! -d "$SDF_OUT" ]; then
   $RTG format -o $SDF_OUT/GCF_016920845.1_GAculeatus_UGA_version5_genomic_shortnames_SDF $REF
fi

###################################################################
#### Mendelian evaluation ####
##############################

MEND_OUTS=$WD/Mendelian_evals/GATK

if [ ! -d "$MEND_OUTS" ]; then
   mkdir -p $MEND_OUTS
fi

PEDIGREE_VCF=/storage/research/iee_evol/DanJ/Stickleback/G_aculeatus/FITNESS/DV_training/Unfiltered_VCF/${CROSS}_pedigree_filtered_GQ20.vcf.gz
PED_FILE=/storage/homefs/dj20y461/Stickleback/G_aculeatus/FITNESS/code/DV_training/PED_files/${CROSS}.ped

$RTG mendelian  -i $PEDIGREE_VCF \
		-t $SDF_OUT \
		-o $MEND_OUTS/${CROSS}.mend_eval_QUAL30.vcf \
		--pedigree $PED_FILE



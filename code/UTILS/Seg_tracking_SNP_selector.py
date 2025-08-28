"""
This script is used to select SNPs based on whether they pass filter thresholds for:
1. Missingness
2. Minor Allele Frequency (MAF)
3. Hardy-Weinberg Equilibrium (HWE)

It outputs a list of SNPs that pass all filters. 

"""

from matplotlib import pyplot as plt
from collections import Counter
import numpy as np
import os
import sys

# define a function to read in the missingness data

def read_lmiss(lmiss_path):
    F_missing = {}

    with open(lmiss_path) as lmiss_file:
        headers = next(lmiss_file)

        for line in lmiss_file:
            CHR, POS, N_DATA, N_GENOTYPE_FILTERED, N_MISS, F_MISS = line.split()[:6]
            F_MISS = float(F_MISS)
            F_missing[(CHR,POS)] = F_MISS

    return F_missing

## define a function to read in the frequencies for each population

def read_freqs(freq_path):
    freqs = {} 
    freqs

    with open(freq_path) as freqs_file:
        headers = next(freqs_file)
        for line in freqs_file:
            chrom, pos, n_alleles, n_chr, freq_1, freq_2 = line.split()[:6]
            allele_1 = freq_1.split(":")[0]
            freq_1 = float(freq_1.split(":")[1])
            allele_2 = freq_2.split(":")[0]
            freq_2 = float(freq_2.split(":")[1])

            freqs[(chrom, pos)] = {}
            freqs[(chrom, pos)]["frequencies"] = (freq_1, freq_2)
            freqs[(chrom, pos)]["alleles"] = (allele_1, allele_2)
            
    return freqs

## define a function to read in the HWE results

def read_HWE(HWE_path):
    
    HWEs = {}

    with open(HWE_path) as HWE_file:
        headers = next(HWE_file)
        for line in HWE_file:
            chrom, pos, obs, exp, ChiSq, P_HWE, P_HET_DEF, P_HET_EXCESS = line.split()[:8]
            HWEs[(chrom, pos)] = {}
            HWEs[(chrom, pos)]["P_HWE"] = float(P_HWE)
            HWEs[(chrom, pos)]["ChiSq"] = float(ChiSq)
            HWEs[(chrom, pos)]["P_HET_DEF"] = float(P_HET_DEF)
            HWEs[(chrom, pos)]["P_HET_EXCESS"] = float(P_HET_EXCESS)
            
    return HWEs


## Get the input arguments
data_dir = sys.argv[1]
suffix = sys.argv[2]
identifier = sys.argv[3] ## eg chromosome, population name etc
groups = sys.argv[4] ## The file identifier to be used to group files for comparison. E.g. population. Comma separated list with no spaces.
segregation_maf_threshold = float(sys.argv[5]) ## The minimum MAF threshold required for a SNP to PASS in a given population

print(groups)
print("searching for files with suffix '%s' and identifier '%s' in directory: '%s'" %(suffix, identifier, data_dir))
print("Grouping files by groups: %s" % groups)

groups = groups.split(",")

## make dictionary to hold file paths for each group
group_file_paths_dict = {group: [] for group in groups}

## get the files in the data directory
for file in os.listdir(data_dir):
    if file.endswith(suffix) and identifier in file:
        #print(file)
        for group in groups:
            if group in file:
                group_file_paths_dict[group].append(os.path.join(data_dir, file))

for group in group_file_paths_dict:
    print("Group '%s' has files: %s" % (group, ", ".join(group_file_paths_dict[group])))

group_freqs_dict = {group: {} for group in groups}

for group in groups:
    print("Reading frequencies for group '%s'" % group)

    for file_path in group_file_paths_dict[group]:
        print("  Reading file: %s" % file_path)
        freqs = read_freqs(file_path)
        group_freqs_dict[group].update(freqs)

print("Finished reading in frequencies for all groups.")

## Assign PASS or FAIL for minimum MAF threshold for all loci, across all groups.

print("Finding SNPs with MAF > %s in all groups" % suffix)

for group in group_freqs_dict:
    for locus in group_freqs_dict[group]:
        maf = min(group_freqs_dict[group][locus]["frequencies"])
        
        if maf >= segregation_maf_threshold:
            group_freqs_dict[group][locus]["MAF_status"] = "PASS"
        else:
            group_freqs_dict[group][locus]["MAF_status"] = "FAIL"

print("Finished assigning PASS/FAIL for MAF threshold for all loci in all groups.")

## Output a per-population summary, and also collect info for the population-specific allele check

cross_group_dict = {}

for group in groups:
    
    with open("%s/%s_SNPs_MAF_PASSFAIL_summary.txt" % (data_dir, group), "w") as out_file:
        
        for locus in group_freqs_dict[group]:

            #print(group, locus)

            ## write summary file

            if all([len(group_freqs_dict[group][locus]["alleles"][0]) == 1,
                    len(group_freqs_dict[group][locus]["alleles"][1]) == 1]):
                out_file.write("%s\t%s\t%s\t%s\t%s\n" % (group, 
                                                         "\t".join([str(i) for i in locus]),
                                                         "\t".join([str(i) for i in group_freqs_dict[group][locus]["frequencies"]]),
                                                         "\t".join([str(i) for i in group_freqs_dict[group][locus]["alleles"]]),
                                                         group_freqs_dict[group][locus]["MAF_status"]))

            if not locus in cross_group_dict:
                cross_group_dict[locus] = {}
            if locus not in group_freqs_dict[group]:
                cross_group_dict[locus][group] = {}
                cross_group_dict[locus][group]["frequencies"] = (0,0)
                cross_group_dict[locus][group]["alleles"] = ("NA","NA")
                cross_group_dict[locus][group]["MAF_status"] = "FAIL"
            else:
                cross_group_dict[locus][group] = group_freqs_dict[group][locus].copy()


    print("%s summary written to %s/%s_SNPs_MAF_PASSFAIL_summary.txt" % (group, data_dir, group))
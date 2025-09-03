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
identifier = sys.argv[2] ## eg chromosome, population name etc
groups = sys.argv[3] ## The file identifier to be used to group files for comparison. E.g. population. Comma separated list with no spaces.
segregation_maf_threshold = float(sys.argv[4]) ## The minimum MAF threshold required for a SNP to PASS in a given population

print(groups)
print("searching for files with identifier '%s' in directory: '%s'" %(identifier, data_dir))
print("Grouping files by groups: %s" % groups)

groups = groups.split(",")

## make dictionary to hold file paths for each group
group_file_paths_dict = {group: [] for group in groups}

## get the files in the data directory
for file in os.listdir(data_dir):
    if all([identifier in file,
            file.endswith(".frq")]):
        
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

print("Finding SNPs with MAF > %s in all groups" % segregation_maf_threshold)

for group in group_freqs_dict:
    for locus in group_freqs_dict[group]:
        maf = min(group_freqs_dict[group][locus]["frequencies"])
        
        if maf >= segregation_maf_threshold:
            group_freqs_dict[group][locus]["MAF_status"] = "PASS"
        else:
            group_freqs_dict[group][locus]["MAF_status"] = "FAIL"

print("Finished assigning PASS/FAIL for MAF threshold for all loci in all groups.")

## Make a locus-wise dictionary to hold info on how many groups a locus passes the MAF threshold in

locus_group_MAF_status = {}

for group in group_freqs_dict:
    for locus in group_freqs_dict[group]:
        if locus not in locus_group_MAF_status:
            locus_group_MAF_status[locus] = {"groups_passed": [], "total_groups": [], "frequencies": []}
        
        locus_group_MAF_status[locus]["total_groups"].append(group)
        locus_group_MAF_status[locus]["frequencies"].append(group_freqs_dict[group][locus]["frequencies"])
        
        if group_freqs_dict[group][locus]["MAF_status"] == "PASS":
            locus_group_MAF_status[locus]["groups_passed"].append(group)

with open(os.path.join(data_dir, "%s_Segregation_SNPs_MAF_%s_PASS_FAIL.txt" % (identifier, segregation_maf_threshold)), "w") as out_file:
    out_file.write("CHROM\tPOS\tN_GROUPS_TOTAL\tALL_GROUPS\tN_GROUPS_PASSED\tGROUPS_PASSED\tALL_FREQUENCIES\n")
    
    for locus in locus_group_MAF_status:
        chrom, pos = locus
        n_groups_total = len(locus_group_MAF_status[locus]["total_groups"])
        n_groups_passed = len(locus_group_MAF_status[locus]["groups_passed"])
        groups_passed = ",".join(locus_group_MAF_status[locus]["groups_passed"])
        total_groups = ",".join(locus_group_MAF_status[locus]["total_groups"])

        #if len(locus_group_MAF_status[locus]["frequencies"]) < 2:
        #print(locus, locus_group_MAF_status[locus]["frequencies"])

        freq_strings = []
        for i in locus_group_MAF_status[locus]["frequencies"]:
            freq_str = "%s,%s" % (i[0], i[1])
            freq_strings.append(freq_str)

        all_freqs = ";".join(freq_strings)
                
        out_file.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (chrom, pos, n_groups_total, total_groups, n_groups_passed, groups_passed, all_freqs))

print("Wrote MAF PASS/FAIL results to file: %s" % os.path.join(data_dir, "%s_Segregation_SNPs_MAF_%s_PASS_FAIL.txt" % (identifier, segregation_maf_threshold)))
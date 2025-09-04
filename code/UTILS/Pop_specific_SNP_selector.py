"""
This script is used to select population-specific SNPs based on frequencies:

It outputs a list of population specific SNPs, their alleles, and frequencies in each population.

"""

from matplotlib import pyplot as plt
from collections import Counter
import numpy as np
import os
import sys

## Get the input arguments
data_dir = sys.argv[1]
identifier = sys.argv[2] ## eg chromosome, population name etc
groups = sys.argv[3] ## The file identifier to be used to group files for comparison. E.g. population. Comma separated list with no spaces.
pop_specific_threshold = float(sys.argv[4]) ## The minimum allele frequency required for a SNP to be considered present in a population to which it is specific 


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


print(groups)
print("searching for files with identifier '%s' in directory: '%s'" %(identifier, data_dir))
print("Grouping files by groups: %s" % groups)

groups = groups.split(",")

## make dictionary to hold file paths for each group
group_file_paths_dict = {group: [] for group in groups}

## get the files in the data directory
for file in os.listdir(data_dir):
    if identifier in file:
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


print("Finding loci with population specific alleles...")

cross_group_dict = {}

for group in groups:
    
    for locus in group_freqs_dict[group]:
            
        if all([len(group_freqs_dict[group][locus]["alleles"][0]) == 1,
                len(group_freqs_dict[group][locus]["alleles"][1]) == 1]):
                
            if not locus in cross_group_dict:
                cross_group_dict[locus] = {}
            cross_group_dict[locus][group] = group_freqs_dict[group][locus]


pop_specific_allele_loci = []

for group in groups:

    others = [i for i in groups if i != group]

    for locus in cross_group_dict:

        found_in_others = False

        if group in cross_group_dict[locus]:

            if cross_group_dict[locus][group]["frequencies"][1] >= pop_specific_threshold:
            
                for other in others:
                    if other in cross_group_dict[locus]:
                        if cross_group_dict[locus][other]["frequencies"][1] > 0: 
                            found_in_others = True   

                if found_in_others == False:
                    pop_specific_allele_loci.append(locus)

print("Number of loci with population specific alleles: %s" % len(pop_specific_allele_loci))

if len(pop_specific_allele_loci) > 0:
    
    print("Writing them to file")

    with open("%s/%s_population_specific_allele_loci.txt" % (data_dir, identifier), "w") as out_file:

        out_file.write("#CHROM\tPOS\t%s\n" % "\t".join(groups))

        for locus in pop_specific_allele_loci:
            output_line = [locus[0], locus[1]]
            output_line.append(",".join(group_freqs_dict[group][locus]["alleles"])) ## hacky, but doesn't matter what group. 
            
            for group in groups:
                #print(group, cross_group_dict[locus])
                if group in cross_group_dict[locus]:
                    output_line.append(cross_group_dict[locus][group]["frequencies"][1])
                else:
                    output_line.append("NA")
            out_file.write("\t".join([str(i) for i in output_line]) + "\n")   

    print("%s/%s_population_specific_allele_loci.txt" % (data_dir, identifier))

else:
    print("No population specific allele loci to output.")

def find_files(indir, suffix):
    
    """
    This function will find all files in a folder with a given suffix.

    """

    import os  
    
    file_list = []

    for file in os.listdir(indir):
        if file.endswith(suffix):
            file_list.append("%s/%s" % (indir,file))
            
    return file_list

### Ok, so for a single file, we need to get autosomal coverage, and then compare to the sex chromosome. 

## list of autosomes (excluding unplaced scaffs)

autos = ["NC_053212.1_chromosome_1","NC_053213.1_chromosome_2","NC_053214.1_chromosome_3","NC_053215.1_chromosome_4","NC_053216.1_chromosome_5",
         "NC_053217.1_chromosome_6","NC_053218.1_chromosome_7","NC_053219.1_chromosome_8","NC_053220.1_chromosome_9","NC_053221.1_chromosome_10",
         "NC_053222.1_chromosome_11","NC_053223.1_chromosome_12","NC_053224.1_chromosome_13","NC_053225.1_chromosome_14","NC_053226.1_chromosome_15",
         "NC_053227.1_chromosome_16","NC_053228.1_chromosome_17","NC_053229.1_chromosome_18","NC_053231.1_chromosome_20",
         "NC_053232.1_chromosome_21"]

sex_chrom = ["NC_053230.1_chromosome_19"]

def compare_covs(covfile, 
                 autosomes = autos,
                 sex_chromosomes = sex_chrom,
                 female_lower = 0.9, 
                 female_upper = 1.1, 
                 male_lower = 0.5,
                 male_upper = 0.6):
    
    """
    This function takes a single coverage file, and divides the average X coverage by the average autosomal coverage. 
    It outputs a file with all data, as well as an automatic classification of XX, XY, or ambiguous. 
    It also produces a pdf with a histogram of X/Autosome coverage ratios for the files processed.
    
    """
 
    import numpy as np
    
    autocovs = []
    sexchromcovs = []
    
    sample_name = covfile.rpartition("/")[2].split(".")[0]
        
    with open(covfile) as covhandle:

        headers = next(covhandle)

        for line in covhandle:

            if line.startswith("NC"):

                chrom = line.split()[0]
                cov = float(line.split()[6])

                if chrom in autosomes:
#                print("YES")
                    autocovs.append(cov)
                
                elif chrom in sex_chromosomes:
                    sexchromcovs.append(cov)

#            else:
#                print("%s not used" % chrom)

        ## get averages

        autoavg = np.round(np.mean(autocovs),2)
        autostdv = np.round(np.std(autocovs),2)
        sexchromavg = np.round(np.mean(sexchromcovs),2)

        ## check sex chrom against autosomes



        ratio = np.round(sexchromavg/autoavg,2)

        if female_lower <= ratio <= female_upper:
            result = "XX"

        elif male_lower <= ratio < male_upper:
            result = "XY" 

        else: 
            result = "ambiguous"
        
    return autoavg, autostdv, sexchromavg, ratio, result


def plot_cov_ratios(cov_dict, 
                    wd,
                    female_lower = 0.9, 
                    female_upper = 1.1, 
                    male_lower = 0.5,
                    male_upper = 0.6):

    from matplotlib import pyplot as plt

    ratios = []

    for sample in cov_dict:

        ratios.append(cov_dict[sample]["ratio"])

    plt.figure(figsize = (10,5))
    hist = plt.hist(ratios, bins = 100, color = "black")

    ## 2n cutoffs
    plt.hlines(0-(max(hist[0]*0.9)*0.02), female_lower, female_upper, linestyle = "-", color = "darkred")
    plt.text(female_lower + ((female_upper-female_lower)/2), 0-(max(hist[0]*0.9)*0.2), "Female range", ha = "center", color = "darkred")

    ## 1n cutoffs
    plt.hlines(0-(max(hist[0]*0.9)*0.02), male_lower, male_upper, linestyle = "-", color = "darkblue")
    plt.text(male_lower + ((male_upper-male_lower)/2), 0-(max(hist[0]*0.9)*0.2), "Male range", ha = "center", color = "darkblue")

    plt.title("X / Autosomal coverage per sample (%s samples)" % len(cov_dict))
    plt.ylabel("N samples")

    plt.savefig("%s/X_vs_Auto_coverage_ratios.pdf" % wd)


##########################
####### RUNNING  #########
##########################

import sys

### Thresholds for ratio of coverage for sex chromsoome to autosomes

FEM_lower = 0.9
FEM_upper = 1.1
MALE_lower = 0.7
MALE_upper = 0.8

## find coverage files
covfiles = find_files(sys.argv[1], sys.argv[2])

cov_dict = {}

with open("%s/X_vs_Auto_coverage_ratios.tsv" % sys.argv[1], 'w') as results_handle:

    results_handle.write("SAMPLE\tmean_AUTO\tstdv_AUTO\tmean_SEXCHR\tRATIO\tRESULT\n")

    for file in covfiles:
    
        sample_name = file.rpartition("/")[2].split(".")[0]
        print(sample_name)

        cov_dict[sample_name] = {}
    
        cov_dict[sample_name]["autoavg"], cov_dict[sample_name]["autostdv"], cov_dict[sample_name]["sexchromavg"], cov_dict[sample_name]["ratio"], cov_dict[sample_name]["result"] = compare_covs(file, female_lower = FEM_lower,female_upper = FEM_upper, male_lower = MALE_lower,male_upper = MALE_upper)

        results_handle.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (sample_name, cov_dict[sample_name]["autoavg"], cov_dict[sample_name]["autostdv"], cov_dict[sample_name]["sexchromavg"], cov_dict[sample_name]["ratio"], cov_dict[sample_name]["result"]))

plot_cov_ratios(cov_dict,
                sys.argv[1],
                female_lower = FEM_lower,
                female_upper = FEM_upper, 
                male_lower = MALE_lower,
                male_upper = MALE_upper)



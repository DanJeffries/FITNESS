#### Plotting results of Sequencing tests for the FITNESS project #####

## Bern ## 

## Trimming, aligning, duplicate removal and final coverage

install.packages("ggplot2")
install.packages("reshape2")

library(ggplot2)
library(reshape2)

#### Raw reads vs concentration

N_reads <- read.delim("/Users/dj20y461/Data_temp/Stickleback/FITNESS/concentrations/Bern/raw_reads.txt", header = T)
concs <- read.delim("/Users/dj20y461/Data_temp/Stickleback/FITNESS/concentrations/Bern/concentrations.txt", header = T)

plotdata <- data.frame(N_reads$N_raw_read_pairs, concs$conc_ng_ul)

plot_data_subset <- plotdata[plotdata$concs.conc_ng_ul < 400,] ## removed two outliers with huge concentrations!

plot(plot_data_subset$N_reads.N_raw_read_pairs, plot_data_subset$concs.conc_ng_ul,
     xlab = "N raw reads",
     ylab = "Conc (ng/ul)",
     main = "Bern raw reads vs. concentration")
abline(lm(concs.conc_ng_ul~N_reads.N_raw_read_pairs,data=plot_data_subset),col='red')


#########################################
## Raw reads vs trimming vs duplicates ##  
########################################


bern_data <- read.delim("~/Downloads/FITNESS_seq_test_Bern_Overview.csv", header = T, sep = ',')

## subset
bern_subset <- data.frame(bern_data$AASample,         ## sample names
                          bern_data$TRIM_dropped,     ## N reads trimmed
                          bern_data$ALIGN_unaligned,  ## N reads that didn't align
                          bern_data$DUPLICATES_total, ## N duplicate reads
                          bern_data$Total_remaining)     ## N reads reamining for SNP calling

## reshape 
bern_sub_reshaped <- melt(bern_subset)

ggplot(bern_sub_reshaped, aes(fill=variable, y=value, x=bern_data.AASample)) + 
  geom_bar(position="stack", stat="identity")

bern_input = mean(bern_data$N_raw_read_pairs, na.rm = T)
bern_remaining_mean = mean(bern_data$Total_remaining, na.rm = T)


### MCGILL ###

#### Raw reads vs concentration

MG_N_reads <- read.delim("/Users/dj20y461/Data_temp/Stickleback/FITNESS/concentrations/Mcgill/raw_reads.txt", header = T)
MG_concs <- read.delim("/Users/dj20y461/Data_temp/Stickleback/FITNESS/concentrations/Mcgill/concentrations.txt", header = T)


MG_plotdata <- data.frame(MG_N_reads$N_raw_read_pairs, MG_concs$conc_ng_µ.)

#plot_data_subset <- plotdata[plotdata$concs.conc_ng_ul < 400,]

plot(MG_plotdata$MG_N_reads.N_raw_read_pairs, MG_plotdata$MG_concs.conc_ng_µ.,
     xlab = "N raw reads",
     ylab = "Conc (ng/ul)",
     main = "McGill raw reads vs. concentration")

abline(lm(MG_concs.conc_ng_µ.~MG_N_reads.N_raw_read_pairs,data=MG_plotdata),col='red')


#########################################
## Raw reads vs trimming vs duplicates ##  
########################################

mcgill_data <- read.delim("~/Downloads/FITNESS_seq_test_McGill_Overview.csv", header = T, sep = ',')

## subset
mcgill_subset <- data.frame(Samples = mcgill_data$AASample,         ## sample names
                          Trimmed = mcgill_data$TRIM_dropped,     ## N reads trimmed
                          Unaligned = mcgill_data$ALIGN_unmapped,  ## N reads that didn't align
                          Duplicates = mcgill_data$DUPLICATES_total, ## N duplicate reads
                          Reads_remaining = mcgill_data$Total_remaining)     ## N reads reamining for SNP calling

## reshape 
mcgill_sub_reshaped <- melt(mcgill_subset)

ggplot(mcgill_sub_reshaped, aes(fill=variable, y=value, x=Samples)) + 
  geom_bar(position="stack", stat="identity")

mcgill_input = mean(mcgill_data$N_raw_read_pairs, na.rm = T)
remaining_mean = mean(mcgill_data$Total_remaining)


##############################################
####### Library complexity vs Conc ###########
##############################################

## Bern

bern_lib_comp <- read.delim("~/Dropbox/Sexy_sticklebacks/FITNESS/TEST_JAN_2023/lib_complexity/Bern_lib_comp_vs_conc.txt", header = T)

plot(bern_lib_comp$Concentration, bern_lib_comp$Estimated_lib_size,
     xlab = "Concentration (ng/ul",
     ylab = "Est. library size (N unique fragments)",
     main = "Bern")

abline(lm(bern_lib_comp$Estimated_lib_size~bern_lib_comp$Concentration),col='red')


## McGill

McGill_lib_comp <- read.delim("~/Dropbox/Sexy_sticklebacks/FITNESS/TEST_JAN_2023/lib_complexity/McGill_lib_comp_vs_conc.txt", header = T)

plot(McGill_lib_comp$Concentration, McGill_lib_comp$Estimated_lib_size,
     xlab = "Concentration (ng/ul",
     ylab = "Est. library size (N unique fragments)",
     main = "McGill")

abline(lm(McGill_lib_comp$Estimated_lib_size~McGill_lib_comp$Concentration),col='red')







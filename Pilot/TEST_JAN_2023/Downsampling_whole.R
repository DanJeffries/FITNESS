##### Comparing average coverage with N SNPs called per population. 

### Coverage

## I subsampled each bam to 60%, 70%, 80%, 90% of the original reads. 
## Thus the coverage, as a fraction of the whole dataset
## should more or less correspond to these proportions . . . 

setwd("~/Data_temp/Stickleback/FITNESS/MCGILL_COV_TESTS/whole/covs/")

### So I will calculate average per pop and per downsampling fraction

## 0.6 
CN06 <- read.delim("./0.6/CN_covs.txt", header = F)
FG06 <- read.delim("./0.6/FG_covs.txt", header = F)
LG06 <- read.delim("./0.6/LG_covs.txt", header = F)
SL06 <- read.delim("./0.6/SL_covs.txt", header = F)
TN06 <- read.delim("./0.6/TN_covs.txt", header = F)
WB06 <- read.delim("./0.6/WB_covs.txt", header = F)
WK06 <- read.delim("./0.6/WK_covs.txt", header = F)
WT06 <- read.delim("./0.6/WT_covs.txt", header = F)

Avg_cov_06 <- c(mean(CN06$V1),
                mean(FG06$V1),
                mean(LG06$V1),
                mean(SL06$V1),
                mean(TN06$V1),
                mean(WB06$V1),
                mean(WK06$V1),
                mean(WT06$V1))

## 0.7 
CN07 <- read.delim("./0.7/CN_covs.txt", header = F)
FG07 <- read.delim("./0.7/FG_covs.txt", header = F)
LG07 <- read.delim("./0.7/LG_covs.txt", header = F)
SL07 <- read.delim("./0.7/SL_covs.txt", header = F)
TN07 <- read.delim("./0.7/TN_covs.txt", header = F)
WB07 <- read.delim("./0.7/WB_covs.txt", header = F)
WK07 <- read.delim("./0.7/WK_covs.txt", header = F)
WT07 <- read.delim("./0.7/WT_covs.txt", header = F)

Avg_cov_07 <- c(mean(CN07$V1),
                mean(FG07$V1),
                mean(LG07$V1),
                mean(SL07$V1),
                mean(TN07$V1),
                mean(WB07$V1),
                mean(WK07$V1),
                mean(WT07$V1))


## 0.8 
CN08 <- read.delim("./0.8/CN_covs.txt", header = F)
FG08 <- read.delim("./0.8/FG_covs.txt", header = F)
LG08 <- read.delim("./0.8/LG_covs.txt", header = F)
SL08 <- read.delim("./0.8/SL_covs.txt", header = F)
TN08 <- read.delim("./0.8/TN_covs.txt", header = F)
WB08 <- read.delim("./0.8/WB_covs.txt", header = F)
WK08 <- read.delim("./0.8/WK_covs.txt", header = F)
WT08 <- read.delim("./0.8/WT_covs.txt", header = F)

Avg_cov_08 <- c(mean(CN08$V1),
                mean(FG08$V1),
                mean(LG08$V1),
                mean(SL08$V1),
                mean(TN08$V1),
                mean(WB08$V1),
                mean(WK08$V1),
                mean(WT08$V1))

## 0.9 
CN09 <- read.delim("./0.9/CN_covs.txt", header = F)
FG09 <- read.delim("./0.9/FG_covs.txt", header = F)
LG09 <- read.delim("./0.9/LG_covs.txt", header = F)
SL09 <- read.delim("./0.9/SL_covs.txt", header = F)
TN09 <- read.delim("./0.9/TN_covs.txt", header = F)
WB09 <- read.delim("./0.9/WB_covs.txt", header = F)
WK09 <- read.delim("./0.9/WK_covs.txt", header = F)
WT09 <- read.delim("./0.9/WT_covs.txt", header = F)

Avg_cov_09 <- c(mean(CN09$V1),
                mean(FG09$V1),
                mean(LG09$V1),
                mean(SL09$V1),
                mean(TN09$V1),
                mean(WB09$V1),
                mean(WK09$V1),
                mean(WT09$V1))

SD_cov_09 <- c(sd(CN09$V1),
                sd(FG09$V1),
                sd(LG09$V1),
                sd(SL09$V1),
                sd(TN09$V1),
                sd(WB09$V1),
                sd(WK09$V1),
                sd(WT09$V1))

## 1.0 
CN1 <- read.delim("./1.0//CN_covs.txt", header = F)
FG1 <- read.delim("./1.0/FG_covs.txt", header = F)
LG1 <- read.delim("./1.0/LG_covs.txt", header = F)
SL1 <- read.delim("./1.0/SL_covs.txt", header = F)
TN1 <- read.delim("./1.0/TN_covs.txt", header = F)
WB1 <- read.delim("./1.0/WB_covs.txt", header = F)
WK1 <- read.delim("./1.0/WK_covs.txt", header = F)
WT1 <- read.delim("./1.0/WT_covs.txt", header = F)

Avg_cov_1 <- c(mean(CN1$V1),
                mean(FG1$V1),
                mean(LG1$V1),
                mean(SL1$V1),
                mean(TN1$V1),
                mean(WB1$V1),
                mean(WK1$V1),
                mean(WT1$V1))


## Additional 0.4

setwd("~/Data_temp/Stickleback/FITNESS/MCGILL_COV_TESTS/Additional/covs/0.4/")

CN_A04 <- read.delim("./CN.covs", sep = " ", header = F)
SL_A04 <- read.delim("./SL.covs", sep = " ", header = F)
TN_A04 <- read.delim("./TN.covs", sep = " ", header = F)
WK_A04 <- read.delim("./WK.covs", sep = " ", header = F)
WT_A04 <- read.delim("./WT.covs", sep = " ", header = F)

Avg_cov_A04 <- c(mean(CN_A04$V1),
                NA,
                NA,
                mean(SL_A04$V1),
                mean(TN_A04$V1),
                NA,
                mean(WK_A04$V1),
                mean(WT_A04$V1))


## Additional 0.6

setwd("~/Data_temp/Stickleback/FITNESS/MCGILL_COV_TESTS/Additional/covs/0.6/")

CN_A06 <- read.delim("./CN.covs", sep = " ", header = F)
SL_A06 <- read.delim("./SL.covs", sep = " ", header = F)
TN_A06 <- read.delim("./TN.covs", sep = " ", header = F)
WK_A06 <- read.delim("./WK.covs", sep = " ", header = F)
WT_A06 <- read.delim("./WT.covs", sep = " ", header = F)

Avg_cov_A06 <- c(mean(CN_A06$V1),
                NA,
                NA,
                mean(SL_A06$V1),
                mean(TN_A06$V1),
                NA,
                mean(WK_A06$V1),
                mean(WT_A06$V1))

## Additional 0.8

setwd("~/Data_temp/Stickleback/FITNESS/MCGILL_COV_TESTS/Additional/covs/0.8/")

CN_A08 <- read.delim("./CN.covs", sep = " ", header = F)
SL_A08 <- read.delim("./SL.covs", sep = " ", header = F)
TN_A08 <- read.delim("./TN.covs", sep = " ", header = F)
WK_A08 <- read.delim("./WK.covs", sep = " ", header = F)
WT_A08 <- read.delim("./WT.covs", sep = " ", header = F)

Avg_cov_A08 <- c(mean(CN_A08$V1),
                NA,
                NA,
                mean(SL_A08$V1),
                mean(TN_A08$V1),
                NA,
                mean(WK_A08$V1),
                mean(WT_A08$V1))

## Additional 1.0

setwd("~/Data_temp/Stickleback/FITNESS/MCGILL_COV_TESTS/Additional/covs/1.0/")

CN_A1 <- read.delim("./CN.covs", sep = " ", header = F)
SL_A1 <- read.delim("./SL.covs", sep = " ", header = F)
TN_A1 <- read.delim("./TN.covs", sep = " ", header = F)
WK_A1 <- read.delim("./WK.covs", sep = " ", header = F)
WT_A1 <- read.delim("./WT.covs", sep = " ", header = F)


Avg_cov_A1 <- c(mean(CN_A1$V2),
                NA,
                NA,
                mean(SL_A1$V2),
                mean(TN_A1$V2),
                NA,
                mean(WK_A1$V2),
                mean(WT_A1$V2))


cov_avgs <- data.frame(O6= Avg_cov_06,
                       O7= Avg_cov_07,
                       O8= Avg_cov_08,
                       O9= Avg_cov_09,
                       O10= Avg_cov_1,
                       A04 = Avg_cov_A04,
                       A06 = Avg_cov_A06,
                       A08 = Avg_cov_A08,
                       A1 = Avg_cov_A1,
                       row.names = c("CN","FG","LG","SL",
                                     "TN","WB","WK","WT"))

## Just checking the coverage improvement . . . 

CN_improve <- mean(CN_A1$V2)/mean(CN1$V1)
SL_improve <- mean(SL_A1$V2)/mean(SL1$V1)
TN_improve <- mean(TN_A1$V2)/mean(TN1$V1)
WK_improve <- mean(WK_A1$V2)/mean(WK1$V1)
WT_improve <- mean(WT_A1$V2)/mean(WT1$V1)

mean(CN_improve,
     SL_improve,
     TN_improve,
     WK_improve,
     WT_improve)

## AVEERAGE OF 3.55 FOLD INCREASE IN COV FROM THE ADDITIONAL SEQUENCING. 


###########################
###### SNP counts #########
###########################

SNPs_CN06 <- 27710
SNPs_FG06 <- 6306
SNPs_LG06 <- 6918
SNPs_SL06 <- 18255
SNPs_TN06 <- 4588
SNPs_WB06 <- 7004
SNPs_WK06 <- 9418
SNPs_WT06 <- 14771

SNPs_CN07 <- 41569
SNPs_FG07 <- 9098
SNPs_LG07 <- 9947
SNPs_SL07 <- 26269
SNPs_TN07 <- 6489
SNPs_WB07 <- 10033
SNPs_WK07 <- 13769
SNPs_WT07 <- 21234

SNPs_CN08 <- 62804
SNPs_FG08 <- 12418
SNPs_LG08 <- 13782
SNPs_SL08 <- 36919
SNPs_TN08 <- 8733
SNPs_WB08 <- 13923
SNPs_WK08 <- 18784
SNPs_WT08 <- 29501

SNPs_CN09 <- 95886
SNPs_FG09 <- 16281
SNPs_LG09 <- 18372
SNPs_SL09 <- 50941
SNPs_TN09 <- 11429
SNPs_WB09 <- 18267
SNPs_WK09 <- 24845
SNPs_WT09 <- 39980

SNPs_CN1 <- 148440
SNPs_FG1 <- 21158
SNPs_LG1 <- 24270
SNPs_SL1 <- 72593
SNPs_TN1 <- 14841
SNPs_WB1 <- 23732
SNPs_WK1 <- 32334
SNPs_WT1 <- 54145

SNPs_CNA04 <- 852840
SNPs_FGA04 <- NA
SNPs_LGA04 <- NA
SNPs_SLA04 <- 2832882
SNPs_TNA04 <- 31242
SNPs_WBA04 <- NA
SNPs_WKA04 <- 345032
SNPs_WTA04 <- 1213279

SNPs_CNA06 <- 4632821
SNPs_FGA06 <- NA
SNPs_LGA06 <- NA
SNPs_SLA06 <- 6779977
SNPs_TNA06 <- 232538
SNPs_WBA06 <- NA
SNPs_WKA06 <- 3591966
SNPs_WTA06 <- 5928895

SNPs_CNA08 <- 6905740
SNPs_FGA08 <- NA
SNPs_LGA08 <- NA
SNPs_SLA08 <- 7516158
SNPs_TNA08 <- 1029243
SNPs_WBA08 <- NA
SNPs_WKA08 <- 6361933
SNPs_WTA08 <- 7183176

SNPs_CNA1 <- 7582121
SNPs_FGA1 <- NA
SNPs_LGA1 <- NA
SNPs_SLA1 <- 7741901
SNPs_TNA1 <- 2605946
SNPs_WBA1 <- NA
SNPs_WKA1 <- 7276214
SNPs_WTA1 <- 7525499


SNPs_all <- data.frame(O6=c(SNPs_CN06,
                            SNPs_FG06,
                            SNPs_LG06,
                            SNPs_SL06,
                            SNPs_TN06,
                            SNPs_WB06,
                            SNPs_WK06,
                            SNPs_WT06),
                       O7=c(SNPs_CN07,
                            SNPs_FG07,
                            SNPs_LG07,
                            SNPs_SL07,
                            SNPs_TN07,
                            SNPs_WB07,
                            SNPs_WK07,
                            SNPs_WT07),
                       O8=c(SNPs_CN08,
                            SNPs_FG08,
                            SNPs_LG08,
                            SNPs_SL08,
                            SNPs_TN08,
                            SNPs_WB08,
                            SNPs_WK08,
                            SNPs_WT08),
                       O9=c(SNPs_CN09,
                            SNPs_FG09,
                            SNPs_LG09,
                            SNPs_SL09,
                            SNPs_TN09,
                            SNPs_WB09,
                            SNPs_WK09,
                            SNPs_WT09),
                       O10=c(SNPs_CN1,
                             SNPs_FG1,
                             SNPs_LG1,
                             SNPs_SL1,
                             SNPs_TN1,
                             SNPs_WB1,
                             SNPs_WK1,
                             SNPs_WT1),
                       A04=c(SNPs_CNA04,
                            SNPs_FGA04,
                            SNPs_LGA04,
                            SNPs_SLA04,
                            SNPs_TNA04,
                            SNPs_WBA04,
                            SNPs_WKA04,
                            SNPs_WTA04),
                       A06=c(SNPs_CNA06,
                            SNPs_FGA06,
                            SNPs_LGA06,
                            SNPs_SLA06,
                            SNPs_TNA06,
                            SNPs_WBA06,
                            SNPs_WKA06,
                            SNPs_WTA06),
                       A08=c(SNPs_CNA08,
                            SNPs_FGA08,
                            SNPs_LGA08,
                            SNPs_SLA08,
                            SNPs_TNA08,
                            SNPs_WBA08,
                            SNPs_WKA08,
                            SNPs_WTA08),
                       A1=c(SNPs_CNA1,
                            SNPs_FGA1,
                            SNPs_LGA1,
                            SNPs_SLA1,
                            SNPs_TNA1,
                            SNPs_WBA1,
                            SNPs_WKA1,
                            SNPs_WTA1),
                        row.names = c("CN","FG","LG","SL",
                                      "TN","WB","WK","WT"))


## And now checking the SNP discovery improvement . . . 

CN_SNP_improve <- SNPs_CNA1/SNPs_CN1
SL_SNP_improve <- SNPs_SLA1/SNPs_SL1
TN_SNP_improve <- SNPs_TNA1/SNPs_TN1
WK_SNP_improve <- SNPs_WKA1/SNPs_WK1
WT_SNP_improve <- SNPs_WTA1/SNPs_WT1

mean(c(CN_SNP_improve,
     SL_SNP_improve,
     TN_SNP_improve,
     WK_SNP_improve,
     WT_SNP_improve))

## AVERAGE OF 3.55 FOLD INCREASE IN COV FROM THE ADDITIONAL SEQUENCING. 


cov_avgs_tr <- t(cov_avgs)
SNPs_all_tr <- t(SNPs_all)

x_min <- min(cov_avgs_tr[c(1,2,3,4,5),], na.rm = T)-(min(cov_avgs_tr[c(1,2,3,4,5),], na.rm = T)/10)
y_min <- min(SNPs_all_tr[c(1,2,3,4,5),], na.rm = T)-(min(SNPs_all_tr[c(1,2,3,4,5),], na.rm = T)/10)

x_max <- max(cov_avgs_tr[c(1,2,3,4,5),], na.rm = T)+(max(cov_avgs_tr[c(1,2,3,4,5),], na.rm = T)/10)
y_max <- max(SNPs_all_tr[c(1,2,3,4,5),], na.rm = T)+(max(SNPs_all_tr[c(1,2,3,4,5),], na.rm = T)/10)


####################################
### PLOTTING previous data tests ###
####################################


plot(cov_avgs_tr[c(1,2,3,4,5),1], SNPs_all_tr[c(1,2,3,4,5),1], 
     ylim = c(y_min,y_max), xlim=c(x_min,x_max), 
     pch=16,
     xlab = "Depth",
     ylab = "N SNPs per population (reads)")
lines(cov_avgs_tr[c(1,2,3,4,5),1], SNPs_all_tr[c(1,2,3,4,5),1])

points(cov_avgs_tr[c(1,2,3,4,5),2], SNPs_all_tr[c(1,2,3,4,5),2], col = "blue", pch=16)
lines(cov_avgs_tr[c(1,2,3,4,5),2], SNPs_all_tr[c(1,2,3,4,5),2], col = "blue")

points(cov_avgs_tr[c(1,2,3,4,5),3], SNPs_all_tr[c(1,2,3,4,5),3], col = "green", pch=16)
lines(cov_avgs_tr[c(1,2,3,4,5),3], SNPs_all_tr[c(1,2,3,4,5),3], col = "green")

points(cov_avgs_tr[c(1,2,3,4,5),4], SNPs_all_tr[c(1,2,3,4,5),4], col = "darkred", pch=16)
lines(cov_avgs_tr[c(1,2,3,4,5),4], SNPs_all_tr[c(1,2,3,4,5),4], col = "darkred")

points(cov_avgs_tr[c(1,2,3,4,5),5], SNPs_all_tr[c(1,2,3,4,5),5], col = "purple", pch=16)
lines(cov_avgs_tr[c(1,2,3,4,5),5], SNPs_all_tr[c(1,2,3,4,5),5], col = "purple")

points(cov_avgs_tr[c(1,2,3,4,5),6], SNPs_all_tr[c(1,2,3,4,5),6], col = "turquoise", pch=16)
lines(cov_avgs_tr[c(1,2,3,4,5),6], SNPs_all_tr[c(1,2,3,4,5),6], col = "turquoise")

points(cov_avgs_tr[c(1,2,3,4,5),7], SNPs_all_tr[c(1,2,3,4,5),7], col = "orange", pch=16)
lines(cov_avgs_tr[c(1,2,3,4,5),7], SNPs_all_tr[c(1,2,3,4,5),7], col = "orange")

points(cov_avgs_tr[,8], SNPs_all_tr[,8], col = "darkgreen", pch=16)
lines(cov_avgs_tr[,8], SNPs_all_tr[,8], col = "darkgreen")

legend(10, y=150000,
       c("CN (N=6)","FG (N=25)","LG (N=28)","SL (N=27)",
         "TN (N=16)","WB (N=29)","WK (N=25)","WT (N=26)"),
       col = c("black","blue", "green", "darkred", "purple", "turquoise", "orange", "darkgreen"),
       pch=16, cex=0.8)




######################################
### PLOTTING additional data tests ###
######################################

x_min <- min(cov_avgs_tr, na.rm = T)-(min(cov_avgs_tr, na.rm = T)/10)
y_min <- min(SNPs_all_tr, na.rm = T)-(min(SNPs_all_tr, na.rm = T)/10)

x_max <- max(cov_avgs_tr, na.rm = T)+(max(cov_avgs_tr, na.rm = T)/10)
y_max <- max(SNPs_all_tr, na.rm = T)+(max(SNPs_all_tr, na.rm = T)/10)


plot(cov_avgs_tr[,1], SNPs_all_tr[,1], 
     ylim = c(y_min,y_max), xlim=c(x_min,x_max), 
     pch=16,
     xlab = "Depth",
     ylab = "N SNPs per population (reads)")
lines(cov_avgs_tr[,1], SNPs_all_tr[,1])

points(cov_avgs_tr[,4], SNPs_all_tr[,4], col = "darkred", pch=16)
lines(cov_avgs_tr[,4], SNPs_all_tr[,4], col = "darkred")

points(cov_avgs_tr[,5], SNPs_all_tr[,5], col = "purple", pch=16)
lines(cov_avgs_tr[,5], SNPs_all_tr[,5], col = "purple")

points(cov_avgs_tr[,7], SNPs_all_tr[,7], col = "orange", pch=16)
lines(cov_avgs_tr[,7], SNPs_all_tr[,7], col = "orange")

points(cov_avgs_tr[,8], SNPs_all_tr[,8], col = "darkgreen", pch=16)
lines(cov_avgs_tr[,8], SNPs_all_tr[,8], col = "darkgreen")

legend(35, y=3000000,
       c("CN (N=6)","FG (N=25)","LG (N=28)","SL (N=27)",
         "TN (N=16)","WB (N=29)","WK (N=25)","WT (N=26)"),
       col = c("black","blue", "green", "darkred", "purple", "turquoise", "orange", "darkgreen"),
       pch=16, cex=0.8)



###### Calculating inbreeding coefficient for each pop to better interpret results

## Doing this just from the total dataset for now.

setwd("~/Data_temp/Stickleback/FITNESS/MCGILL_COV_TESTS/FIS/")

read.delim("1.0/CN_inbreeding.het")










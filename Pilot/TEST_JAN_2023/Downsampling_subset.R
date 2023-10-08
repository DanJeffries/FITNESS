install.packages("vcf2sfs")

## VCFtools filters:

# --remove-indels \
# --max-alleles 2 \
# --max-missing-count 5 \
# --minQ 30 \
# --minGQ 30 \
# 

reads_M <- c(15,19,23,27)
SNPs <- c(4570,8713,14583,23453)

dev.off()
plot(reads_M, SNPs, pch = 16)


# --remove-indels \
# --max-alleles 2 \
# --max-missing-count 5 \
# --minQ 30 \
# --minGQ 30 \
# 

SNPs <- c(13943,25814,43745,75629)

plot(reads_M, SNPs, pch = 16)


###### Plotting SFS #####

source("~/Software/vcf2sfs/vcf2sfs.r")

dev.off()
par(mfrow = c(4,1) , mai = c(0.25, 0.5, 0.15, 0.5)) ## B, L, T, R

x8_GT <- vcf2gt("~/Data_temp/Stickleback/FITNESS/MCGILL_COV_TESTS/x8/McGill_filtered_15.vcf.gz.recode.vcf",
                "~/Data_temp/Stickleback/FITNESS/MCGILL_COV_TESTS/popmap.txt")
sfs_x8 <- gt2sfs.raw(x8_GT, "WT")
plot.sfs(sfs_x8)


x10_GT <- vcf2gt("~/Data_temp/Stickleback/FITNESS/MCGILL_COV_TESTS/x10/McGill_filtered_19.vcf.gz.recode.vcf",
                "~/Data_temp/Stickleback/FITNESS/MCGILL_COV_TESTS/popmap.txt")
sfs_x10 <- gt2sfs.raw(x10_GT, "WT")
plot.sfs(sfs_x10)


x12_GT <- vcf2gt("~/Data_temp/Stickleback/FITNESS/MCGILL_COV_TESTS/x12/McGill_filtered_23.vcf.gz.recode.vcf",
                "~/Data_temp/Stickleback/FITNESS/MCGILL_COV_TESTS/popmap.txt")
sfs_x12 <- gt2sfs.raw(x12_GT, "WT")
plot.sfs(sfs_x12)


x14_GT <- vcf2gt("~/Data_temp/Stickleback/FITNESS/MCGILL_COV_TESTS/x14/McGill_filtered_27.vcf.gz.recode.vcf",
                "~/Data_temp/Stickleback/FITNESS/MCGILL_COV_TESTS/popmap.txt")
sfs_x14 <- gt2sfs.raw(x14_GT, "WT")
plot.sfs(sfs_x14)

#### heterozygosity vs coverage ####

# is there a correlation between het and coverage?#

## x8

covs <- read.delim("~/Data_temp/Stickleback/FITNESS/MCGILL_COV_TESTS/x8/STATS.ldepth.mean")
allhets <- read.delim("~/Data_temp/Stickleback/FITNESS/MCGILL_COV_TESTS/x8/STATS.hwe")
hets <- read.delim("~/Data_temp/Stickleback/FITNESS/MCGILL_COV_TESTS/x8/STATS.obshet")

length(covs$CHROM)
length(hets$CHR)

head(hets)
hets$INDIV <- (hets$OBS.HOM1 + hets$HET + hets$HOM2.)
hets$PROP <- (hets$HET/hets$INDIV)


plotstats <- data.frame(hets$CHR, hets$POS, hets$PROP, allhets$P_HWE, covs$MEAN_DEPTH)
head(plotstats)

str(HWE_stats)

HWE_stats <- subset(plotstats, allhets.P_HWE > 0.05)
HWE_stats
dev.off()

plot(HWE_stats$hets.PROP, HWE_stats$covs.MEAN_DEPTH, ylim = c(0,100))
abline(lm(covs.MEAN_DEPTH~hets.PROP,data=HWE_stats),col='red')


## x10

covs <- read.delim("~/Data_temp/Stickleback/FITNESS/MCGILL_COV_TESTS/x10/STATS.ldepth.mean")
allhets <- read.delim("~/Data_temp/Stickleback/FITNESS/MCGILL_COV_TESTS/x10/STATS.hwe")
hets <- read.delim("~/Data_temp/Stickleback/FITNESS/MCGILL_COV_TESTS/x10/STATS.obshet")

length(covs$CHROM)
length(hets$CHR)

head(hets)
hets$INDIV <- (hets$HOM1 + hets$HET + hets$HOM2)
hets$PROP <- (hets$HET/hets$INDIV)


plotstats <- data.frame(hets$CHR, hets$POS, hets$PROP, allhets$P_HWE, covs$MEAN_DEPTH)
head(plotstats)

str(HWE_stats)

HWE_stats <- subset(plotstats, allhets.P_HWE > 0.05)
HWE_stats
dev.off()

plot(HWE_stats$hets.PROP, HWE_stats$covs.MEAN_DEPTH, ylim = c(0,100))
abline(lm(covs.MEAN_DEPTH~hets.PROP,data=HWE_stats),col='red')

lm(covs.MEAN_DEPTH~hets.PROP,data=HWE_stats)









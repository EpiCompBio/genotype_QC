#############################
# Plink output scripts

# To run after extracting SNP of interest from plink analysis with --recodeAD function.
# See:
# http://pngu.mgh.harvard.edu/~purcell/plink/dataman.shtml#recode
# http://pngu.mgh.harvard.edu/~purcell/plink/tutorial.shtml#t11
#
#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
setwd('/ifs/projects/proj043/analysis.dir/association_analysis_1.dir/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output",".txt", sep=""), open='a')
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))

##TO DO extract parameters:

# Load a previous R session, data and objects:
#load('R_session_saved_image.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
#R_session_saved_image <- paste('R_session_saved_image_normalisation', '.RData', sep='')
R_session_saved_image_full <- paste('R_session_saved_image', '.RData', sep='')

#############################


#############################
## Update packages if necessary and load them:

#source("http://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
#install.packages('')
#detach("package:pryr", unload=TRUE)
#biocLite('SNPRelate')
#biocLite("GeneticTools")

library(SNPRelate)
library(ggplot2)
#library(GeneticTools)
library(plyr)

# SNPRelate webpages 
# https://github.com/zhengxwen/SNPRelate
# http://corearray.sourceforge.net/tutorials/SNPRelate/

# https://www.biostars.org/p/83232/


#############################
# Read in data with SNP of interest:
# Extract SNPs with plink using:
# plink --noweb --bfile P140343-Results_FinalReport_clean-individuals_and_SNPs \
# --snp rs2282679 --recodeAD --out rs2282679_plink_results

'rs1993116_plink_results.raw'
'rs2060793_plink_results.raw'
'rs2282679_plink_results.raw'
'rs3829251_plink_results.raw'

file_snp <- 'rs1993116_plink_results.raw'
SNP <- 'rs1993116'

SNP_of_interest <- read.table(file_snp, header=T, sep=' ')
head(SNP_of_interest)
head(SNP_of_interest$FID)
head(SNP_of_interest[7])

# Recode SNPs:
SNP_of_interest$genotype[SNP_of_interest$rs2282679_G==0 & SNP_of_interest$rs2282679_HET==0] <- 'GG'
SNP_of_interest$genotype[SNP_of_interest$rs2282679_G==1 & SNP_of_interest$rs2282679_HET==1] <- 'GT'
SNP_of_interest$genotype[SNP_of_interest$rs2282679_G==2 & SNP_of_interest$rs2282679_HET==0] <- 'TT'
SNP_of_interest$genotype[SNP_of_interest$rs2282679_G==NA & SNP_of_interest$rs2282679_HET==NA] <- 'NA'

head(SNP_of_interest)
summary(SNP_of_interest)
count(SNP_of_interest$genotype)
#View(SNP_of_interest)

# Read in phenotype data. Replace -99999 (NA's for plink) to NaN for R:
pheno_data <- read.table('BEST-D_alternate_phenotypes_kit_ID_twice_recoded.txt', 
                              header=T, sep=' ', na.strings = '-99999')
head(pheno_data)
head(pheno_data$FID)
#View(pheno_data)

# Merge with phenotype data:
merged_data <- merge(SNP_of_interest, pheno_data, by='FID')
head(merged_data)
#View(merged_data)

# Recode arm where 0=4000 IU, 1=2000 IU, 2=Placebo:
merged_data$arm2[merged_data$arm==0] <- '4000_IU'
merged_data$arm2[merged_data$arm==1] <- '2000_IU'
merged_data$arm2[merged_data$arm==2] <- 'Placebo'
count(merged_data$arm2)

# Boxplot by genotype and phenotype of interest:
png('boxplot_by_genotype_baseline.png', width = 4, height = 4, units = 'in', res = 300)
boxplot(vitd0~genotype, merged_data, main="Levels by genotype",
        xlab="SNP genotype", ylab="25OHD levels")
dev.off()

png('boxplot_by_genotype_VD6.png', width = 4, height = 4, units = 'in', res = 300)
boxplot_genotype <- ggplot(aes(y = vitd6, x = genotype),
                           data = merged_data) + geom_boxplot(
                           ) + geom_jitter(shape=16, position=position_jitter(0.2)
                           ) + scale_color_brewer(
                             palette="Dark2")
boxplot_genotype
dev.off()

png(paste(SNP, '_boxplot_by_genotype_VD12.png', sep=''), width = 4, height = 4, units = 'in', res = 300)
boxplot_genotype <- ggplot(aes(y = vitd12, x = genotype),
                           data = merged_data) + geom_boxplot(
                             ) + geom_jitter(shape=16, position=position_jitter(0.2)
                                             ) + scale_color_brewer(
                                                palette="Dark2")
boxplot_genotype
dev.off()

png('boxplot_by_genotype_arm.png', width = 6, height = 6, units = 'in', res = 300)
Group <- factor(merged_data$arm2, levels=c("Placebo", "2000_IU", "4000_IU"))
Genotype <- factor(merged_data$genotype, levels=c("TT", "GT", "GG"))
boxplot_genotype_arm <- ggplot(aes(y=vitd12, x=Genotype, fill=Group),
                           data = merged_data) + labs (main=SNP, y='25OHD plasma levels (nmol/L)') + geom_boxplot(
                           ) + scale_color_brewer(palette="Dark2")
boxplot_genotype_arm
dev.off()


png(paste(SNP, '_boxplot_by_genotype_arm.png', sep=''), width = 6, height = 6, units = 'in', res = 300)
Group <- factor(merged_data$arm2, levels=c("Placebo", "2000_IU", "4000_IU"))
#Genotype <- factor(merged_data$genotype, levels=c("TT", "GT", "GG"))
Genotype <- factor(merged_data[7]), levels=c("0", "1", "2"))
boxplot_genotype_arm <- ggplot(aes_string(y=vitd12, x=Genotype, fill=Group),
                               data = merged_data) + labs (y='25OHD plasma levels (nmol/L)') + geom_boxplot(
                               ) + scale_color_brewer(palette="Dark2")
boxplot_genotype_arm
dev.off()

#geom_jitter(shape=16, position=position_jitter(0.2))
#geom_dotplot(binaxis='y', stackdir='down', dotsize=0.3)

#############################
# The end:
# Remove objects that are not necessary to save:
#ls()
#object_sizes <- sapply(ls(), function(x) object.size(get(x)))
#as.matrix(rev(sort(object_sizes))[1:10])

#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image_full, compress='gzip')

#objects_to_save <- (c('normalised_expressed', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

q()

# Next: run the script for xxx.
#############################
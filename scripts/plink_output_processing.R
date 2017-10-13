#############################
# Plink output scripts

# To run after extracting SNP of interest from plink analysis with --recodeAD function.
# See:
# http://zzz.bwh.harvard.edu/plink/dataman.shtml
# http://pngu.mgh.harvard.edu/~purcell/plink/tutorial.shtml#t11
# Illumina files:
# https://support.illumina.com/downloads/humanomniexpress-24-v1-0-product-files.html
# https://support.illumina.com/array/array_kits/humanomniexpress-24-beadchip-kit/downloads.html
# https://support.illumina.com/content/dam/illumina-marketing/documents/products/datasheets/datasheet_human_omni_express.pdf
# See dbSNP and UCSC for genotype and allele info:
# https://www.ncbi.nlm.nih.gov/books/NBK21088/
#############################


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd('/ifs/projects/proj043/analysis.dir/association_analysis_1.dir/')
# setwd('/Users/antoniob/Documents/quickstart_projects/BEST_D_molecular.p_q/results/repro_re_runs/')

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
# load('../../data/raw/R_session_saved_image_xxx.RData')

# Filename to save current R session, data and objects at the end:
#R_session_saved_image_full <- paste('R_session_saved_image_xxx', '.RData', sep='')
#############################


#############################
## Update packages if necessary and load them:

#source("http://bioconductor.org/biocLite.R")
#biocLite("BiocUpgrade")
#install.packages('')
#detach("package:pryr", unload=TRUE)
#biocLite('SNPRelate')
#biocLite("GeneticTools")

# library(SNPRelate)
library(ggplot2)
#library(GeneticTools)
library(plyr)
library(dplyr)
library(cowplot)

# SNPRelate webpages 
# https://github.com/zhengxwen/SNPRelate
# http://corearray.sourceforge.net/tutorials/SNPRelate/

# https://www.biostars.org/p/83232/

source('../../code/BEST_D_molecular/utilities/ggtheme.R')
#############################

#############################
# Read in data with SNP of interest:
# Extract SNPs with plink using:
# plink --noweb --bfile P140343-Results_FinalReport_clean-individuals_and_SNPs \
# --snp rs2282679 --recodeAD --out rs2282679_plink_results

# 'rs1993116_plink_results.raw'
# 'rs2060793_plink_results.raw'
# 'rs2282679_plink_results.raw'
# 'rs3829251_plink_results.raw'

file_snp <- '../../data/raw/rs7041_plink_results.raw'
SNP <- 'rs7041'
gene <- 'GC'

SNP_of_interest <- read.table(file_snp, header=T, sep=' ')
head(SNP_of_interest)
head(SNP_of_interest$FID)
head(SNP_of_interest[7])
# View(SNP_of_interest)
#############################

#############################
# Recode SNPs:
# plink --file data --recodeAD gives:
# http://zzz.bwh.harvard.edu/plink/dataman.shtml
# the default for the additive recoding is to count the number of minor alleles per person
# which, assuming C is the minor allele, will recode genotypes as follows:
#   SNP       SNP_A ,  SNP_HET
# ---       -----    -----
#   A A   ->    0   ,   0
# A C   ->    1   ,   1
# C C   ->    2   ,   0
# 0 0   ->   NA   ,  NA
# For rs7041:
# C	is minor in 1000G, ref observed in dbSNP 150 is G/T (REV), C/A (FWD)
# TOP/BOT strands for Illumina:
# https://www.illumina.com/documents/products/technotes/technote_topbot.pdf
# Illumina manisfest file 'HumanOmniExpress-24-v1-0-B.csv' for rs7041:
# Allele1 - Top	= A ; Allele2 - Top = C
# Plink allele_A meaning

minor <- 'C'
major <- 'A'
# Homozygous major:
SNP_of_interest$genotype[SNP_of_interest$rs7041_A==0 & SNP_of_interest$rs7041_HET==0] <- paste(major, major, sep = '')
# Heterozygous:
SNP_of_interest$genotype[SNP_of_interest$rs7041_A==1 & SNP_of_interest$rs7041_HET==1] <- paste(minor, major, sep = '')
# Homozygous minor:
SNP_of_interest$genotype[SNP_of_interest$rs7041_A==2 & SNP_of_interest$rs7041_HET==0] <- paste(minor, minor, sep = '')
# Missing:
SNP_of_interest$genotype[SNP_of_interest$rs7041_A==NA & SNP_of_interest$rs7041_HET==NA] <- 'NA'

head(SNP_of_interest)
nrow(SNP_of_interest)
summary(SNP_of_interest)
plyr::count(SNP_of_interest$genotype)
sum(plyr::count(SNP_of_interest$genotype)$freq)
#View(SNP_of_interest)
#############################

#############################
# Read in phenotype data. Replace -99999 (NA's for plink) to NaN for R:
pheno_data <- read.table('../../data/raw/BEST-D_alternate_phenotypes_kit_ID_twice_recoded_Txsubgroups_new_variables.txt', 
                         header = T, sep = '\t', na.strings = '-99999')
# First plot was with this phenotype file, only new variables and subgroups should have been added:
# pheno_data <- read.table('BEST-D_alternate_phenotypes_kit_ID_twice_recoded.txt', 
#                               header=T, sep=' ', na.strings = '-99999')
head(pheno_data)
head(pheno_data$FID)
#View(pheno_data)

# Merge with phenotype data:
merged_data <- merge(SNP_of_interest, pheno_data, by='FID')
head(merged_data)
merged_data[1:5, c('FID', 'genotype', 'vitd0', 'vitd12')]
summary(merged_data[, c('FID', 'genotype', 'vitd0', 'vitd12')])
# Check means of VD per genotype:
merged_data[, c('genotype', 'vitd0', 'vitd12')] %>% 
  group_by(genotype) %>% 
  summarise_all(mean, na.rm = TRUE)

#View(merged_data)

# Recode arm where 0=4000 IU, 1=2000 IU, 2=Placebo:
merged_data$arm2[merged_data$arm==0] <- '4000_IU'
merged_data$arm2[merged_data$arm==1] <- '2000_IU'
merged_data$arm2[merged_data$arm==2] <- 'Placebo'
plyr::count(merged_data$arm2)
#############################

#############################
# Exploratory:
# Boxplot by genotype and phenotype of interest:
# svg(paste(SNP, '_boxplot_by_genotype_baseline.svg', sep=''))
boxplot(vitd0~genotype, merged_data,
        main = "Levels by genotype",
        xlab = sprintf("SNP genotype - %s", SNP),
        ylab = "25OHD levels at baseline")
# dev.off()

# svg('boxplot_by_genotype_VD6.svg')
ggplot(aes(y = vitd6, x = genotype), data = merged_data) + 
  geom_boxplot() + 
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  scale_color_brewer(palette="Dark2") +
  theme_Publication()
# dev.off()

# svg(paste(SNP, '_boxplot_by_genotype_VD12.svg', sep=''))
ggplot(aes(y = vitd12, x = genotype), data = merged_data) + 
  geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  scale_color_brewer(palette="Dark2") +
  theme_Publication()
# dev.off()
#############################


#############################
# Actual plots for publication, with and without baseline plot x-axis as genotype:
Group <- factor(merged_data$arm2, levels = c("Placebo", "2000_IU", "4000_IU"), 
                labels = c("Placebo", "2000 IU", "4000 IU"))
hom_minor <- paste(minor, minor, sep = '')
hom_major <- paste(major, major, sep = '')
het <- paste(minor, major, sep = '')
Genotype <- factor(merged_data$genotype, levels = c(hom_minor, het, hom_major),
                   labels = c(hom_minor, het, hom_major))
# At baseline:
timepoint <- 'baseline'
p0 <- ggplot(aes(y = vitd0, x = Genotype, fill = Group), data = merged_data) +
  geom_boxplot(position = position_dodge(1)) +
  scale_color_brewer(palette="Dark2") +
  theme_Publication() +
  theme(legend.title = element_blank()) +
  labs (y = sprintf('25OHD plasma levels (nmol/L) at %s', timepoint),
        x = sprintf('%s (%s)', SNP, gene)) +
  theme(legend.position = "none")
# At 6 months:
timepoint <- '6 months'
p1 <- ggplot(aes(y = vitd6, x = Genotype, fill = Group), data = merged_data) +
  geom_boxplot(position = position_dodge(1)) +
  scale_color_brewer(palette="Dark2") +
  theme_Publication() +
  theme(legend.title = element_blank()) +
  labs (y = sprintf('25OHD plasma levels (nmol/L) at %s', timepoint),
        x = sprintf('%s (%s)', SNP, gene)) +
  theme(legend.position = "bottom")
  # geom_jitter(shape = 16, position = position_jitter(0.1))
  # geom_dotplot(binaxis='y', stackdir='down', dotsize=0.3)
# Get legend:
the_legend <- cowplot::get_legend(p1)
# Re-plot without legend:
p1 <- p1 + theme(legend.position = 'none') 
# At 12 months:
timepoint <- '12 months'
p2 <- ggplot(aes(y = vitd12, x = Genotype, fill = Group), data = merged_data) +
  geom_boxplot(position = position_dodge(1)) + 
  scale_color_brewer(palette = "Dark2") +
  theme_Publication() +
  theme(legend.position = 'none') +
  labs (y = sprintf('25OHD plasma levels (nmol/L) at %s', timepoint),
        x = sprintf('%s (%s)', SNP, gene))
# Get rid of x axis titles:
p0 <- p0 + theme(axis.title.x = element_blank())
p1 <- p1 + theme(axis.title.x = element_blank())
p2 <- p2 + theme(axis.title.x = element_blank())
# Without baseline:
# Use cowplot to arrange plots in one and share a legend:
plot_1 <- plot_grid(p1, NULL, p2,
                    nrow = 1,
                    align = 'vh',
                    labels = c("A", "", "B"),
                    rel_widths = c(1, 0.05, 1))
# Single common x-axis:
plot_1_axis <- ggdraw(add_sub(plot_1,
                                 fontface = 'bold',
                                 size = 13,
                                 label = sprintf('%s (%s)', SNP, gene)
                                 ))
# Final:
plot_1_f <- plot_grid(plot_1_axis, 
                         the_legend,
                         ncol = 1,
                         rel_heights = c(1, 0.05))
# Save the actual plot here, use base graphics device, cowplot's defaults don't look nice here:
plot_name <- sprintf('%s_boxplot_arm_6_and_12_months.svg', SNP)
# A4 paper measures 210 × 297 millimeters or 8.27 × 11.69 inches
svg(plot_name, width = 10, height = 5)
plot_1_f
dev.off()

# Add plot with baseline included:
plot_base <- plot_grid(p0, NULL, p1, NULL, p2,
                    nrow = 1,
                    align = 'vh',
                    labels = c("A", "", "B", "", "C"),
                    rel_widths = c(1, 0.05, 1, 0.05, 1))
# Single common x-axis:
plot_base_axis <- ggdraw(add_sub(plot_base,
                              fontface = 'bold',
                              size = 13,
                              label = sprintf('%s (%s)', SNP, gene)
))
# Final:
plot_base_f <- plot_grid(plot_base_axis,
                         the_legend,
                         ncol = 1,
                         rel_heights = c(1, 0.05))
# Save the plot:
plot_name <- sprintf('%s_boxplot_arm_base_6_and_12_months.svg', SNP)
svg(plot_name, width = 10, height = 5)
plot_base_f
dev.off()
# save_plot(sprintf('%s_cowplot.svg', plot_name), plot_f)
#############################

#############################
# Actual plots for publication, with and without baseline plot x-axis as allocation:
# Only switching x and fill with Group and Genotype as defined above
# At baseline:
timepoint <- 'baseline'
p0 <- ggplot(aes(y = vitd0, x = Group, fill = Genotype), data = merged_data) +
  geom_boxplot(position = position_dodge(1)) +
  scale_color_brewer(palette="Dark2") +
  theme_Publication() +
  theme(legend.title = element_blank()) +
  labs (y = sprintf('25OHD plasma levels (nmol/L) at %s', timepoint),
        x = sprintf('%s (%s)', SNP, gene)) +
  theme(legend.position = "none")
# At 6 months:
timepoint <- '6 months'
p1 <- ggplot(aes(y = vitd6, x = Group, fill = Genotype), data = merged_data) +
  geom_boxplot(position = position_dodge(1)) +
  scale_color_brewer(palette="Dark2") +
  theme_Publication() +
  theme(legend.title = element_blank()) +
  labs (y = sprintf('25OHD plasma levels (nmol/L) at %s', timepoint),
        x = sprintf('%s (%s)', SNP, gene)) +
  theme(legend.position = "bottom")
# Get legend:
the_legend <- cowplot::get_legend(p1)
# Re-plot without legend:
p1 <- p1 + theme(legend.position = 'none') 
# At 12 months:
timepoint <- '12 months'
p2 <- ggplot(aes(y = vitd12, x = Group, fill = Genotype), data = merged_data) +
  geom_boxplot(position = position_dodge(1)) + 
  scale_color_brewer(palette = "Dark2") +
  theme_Publication() +
  theme(legend.position = 'none') +
  labs (y = sprintf('25OHD plasma levels (nmol/L) at %s', timepoint),
        x = sprintf('%s (%s)', SNP, gene))
# Get rid of x axis titles:
p0 <- p0 + theme(axis.title.x = element_blank())
p1 <- p1 + theme(axis.title.x = element_blank())
p2 <- p2 + theme(axis.title.x = element_blank())
# Angle the x-axis text:
# p0 <- p0 + theme(axis.text.x = element_text(angle = -10, vjust = 1))
# p1 <- p1 + theme(axis.text.x = element_text(angle = -10, vjust = 1))
# p2 <- p2 + theme(axis.text.x = element_text(angle = -10, vjust = 1))
# Use cowplot to arrange plots in one and share a legend:
# Without baseline:
plot_1 <- plot_grid(p1, NULL, p2,
                    nrow = 1,
                    align = 'vh',
                    labels = c("A", "", "B"),
                    rel_widths = c(1, 0.05, 1))
# Add a single x-axis title:
plot_1_axis <- ggdraw(add_sub(plot_1,
                                 fontface = 'bold',
                                 size = 13,
                                 label = sprintf('%s (%s)', SNP, gene)
                                 # vpadding = grid::unit(0, "lines"),
                                 # y = 6, x = 0.5, vjust = 4.5)
                              ))
plot_1_axis_f <- plot_grid(plot_1_axis,
                           the_legend,
                           ncol = 1,
                           rel_heights = c(1, 0.05))
# Save the actual plot here, use base graphics device, cowplot's defaults don't look nice here:
plot_name <- sprintf('%s_boxplot_x-axis_arm_6_and_12_months.svg', SNP)
svg(plot_name, width = 10, height = 5)
plot_1_axis_f
dev.off()
# Add plot with baseline included:
plot_base <- plot_grid(p0, NULL, p1, NULL, p2,
                       nrow = 1,
                       align = 'vh',
                       labels = c("A", "", "B", "", "C"),
                       rel_widths = c(1, 0.05, 1, 0.05, 1))
plot_base_axis <- ggdraw(add_sub(plot_base,
                              fontface = 'bold',
                              size = 13,
                              label = sprintf('%s (%s)', SNP, gene)
                              # vpadding = grid::unit(0, "lines"),
                              # y = 6, x = 0.5, vjust = 4.5)
))
plot_base_f <- plot_grid(plot_base_axis, 
                         the_legend,
                         ncol = 1,
                         rel_heights = c(1, 0.05))
# Save the plot:
plot_name <- sprintf('%s_boxplot_x-axis_arm_base_6_and_12_months.svg', SNP)
svg(plot_name, width = 10, height = 5)
plot_base_f
dev.off()
#############################


#############################
# The end:
# Remove objects that are not necessary to save:
#ls()
#object_sizes <- sapply(ls(), function(x) object.size(get(x)))
#as.matrix(rev(sort(object_sizes))[1:10])

#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))

# To save R workspace with all objects to use at a later time:
# save.image(file=R_session_saved_image_full, compress='gzip')

#objects_to_save <- (c('normalised_expressed', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')

sessionInfo()
q()

# Next: run the script for xxx.
#############################
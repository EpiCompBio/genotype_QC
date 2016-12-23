#############################################
# Match metabolomics data and vcf columns
# Antonio J Berlanga-Taylor
# 7 June 2016
# Input: 
# Outputs: 
#############################################


#############################################
##Set working directory and file locations and names of required inputs:
options(echo = TRUE)

# Working directory:
# setwd('/Users/antoniob/Desktop/Airwave_inflammation/results_2.dir/')

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.

output_file <- file(paste("R_session_output_trans_in_cis_counts",Sys.Date(),".txt", sep=""))
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))
getwd()

##TO DO extract parameters:

# Re-load a previous R session, data and objects:
#load('R_session_saved_image_order_and_match.RData', verbose=T)
#load('R_session_saved_image_eQTL_responseQTLs.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_vcf_header_match','.RData', sep='')
R_session_saved_image
####################


####################
# Load packages:
library(data.table)
####################

####################
# Run with cmd arguments:
args <- commandArgs(trailingOnly = TRUE)

# File names:
header_file1 <- as.character(args[1])
# header_file1 <- 'vcf_columns.txt'

header_file2 <- as.character(args[2])
# header_file2 <- 'Airwave_CPMG_Plasma.txt_sample_names_FID_IID.txt'

print(args)
####################

####################
# Set up file naming for outputs:
# to create: 2000+4000-baseline-1.trans_in_cis
# eQTL_file1_base <- strsplit(eQTL_file1, '[.]')
# eQTL_file1_base <- eQTL_file1_base[[1]][1]
# eQTL_file1_base
# eQTL_file1_ext <- strsplit(eQTL_file1, '_')
# eQTL_file1_ext <- eQTL_file1_ext[[1]][2]
# eQTL_file1_ext
# 
# eQTL_file2_base <- strsplit(eQTL_file2, '[.]')
# eQTL_file2_base <- eQTL_file2_base[[1]][1]
# eQTL_file2_base
# eQTL_file2_ext <- strsplit(eQTL_file2, '_')
# eQTL_file2_ext <- eQTL_file2_ext[[1]][2]
# eQTL_file2_ext
# 
# output_file_name <- sprintf('%s_in_%s_%s.eQTL', eQTL_file1_ext, eQTL_file2_ext, eQTL_file1_base)
# output_file_name
####################


####################
# Read data:
# header_data1 <- fread(header_file1, sep = '\t', header = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
header_data1 <- read.csv(header_file1, sep = '\t', header = FALSE, stringsAsFactors = FALSE, strip.white = TRUE, skip = 9)
dim(header_data1)
colnames(header_data1)
header_data1[1, 1:5]
# header_data1[1:2, 1:5, with = F]
# header_data1[1:2, 10:15, with = F]

header_data2 <- read.csv(header_file2, sep = '\t', header = FALSE, stringsAsFactors = FALSE, strip.white = TRUE)
dim(header_data2)
colnames(header_data2)
head(header_data2)
# setkey(header_data2, V1)

# tables()
####################


####################
# Transpose vcf header to allow merging:
header_data1 <- transpose(header_data1)
class(header_data1)
head(header_data1)
# setkey(header_data1, V1)
dim(header_data1)

# tables()
# key(header_data1)
# key(header_data2)
# header_data1
# header_data2
str(header_data1)
str(header_data2)

# setkey(header_data1)
# setkey(header_data2)

# Merge columns:
# merged_headers <- header_data1[header_data2, nomatch = 0] #data.table inner join, returns rows with a match in both tables
merged_headers <- merge(header_data1, header_data2) #data.table inner join, returns rows with a match in both tables
head(merged_headers)
dim(merged_headers)
dim(header_data2)[1] - dim(merged_headers)[1] # Missing samples
# 988 samples with metabolomics data are missing genotype data
####################

####################
# Get columns from vcf file that match metabolomics samples:
class(header_data1[, 1])
# vcf_cols_to_extract <- which(header_data1[, V1] %in% header_data2[, V1])
vcf_cols_to_extract <- which(header_data1[, 'V1'] %in% header_data2[, 'V1'])
length(vcf_cols_to_extract)
t(vcf_cols_to_extract)
head(vcf_cols_to_extract)
# header_data1[1:5, 1, with = F]

# Format to pass to awk/cut:
class(vcf_cols_to_extract)
write.csv(t(vcf_cols_to_extract), 'vcf_cols_to_extract.csv', row.names = FALSE, col.names = FALSE)
# Delete headers:
system('head -n 1 vcf_cols_to_extract.csv')
system("sed '1d' vcf_cols_to_extract.csv > vcf_cols_to_extract_2.csv")
system('head vcf_cols_to_extract_2.csv') # Needs columns 1-9 added to get chrom, pos, info, etc.
system('rm -f vcf_cols_to_extract.csv')

# Test on subset of samples:
# See subset_samples.sh
####################

####################
# Get headers after subsetting:
# system('rsync -vvazPui --delete aberlang@login.cx1.hpc.ic.ac.uk:/groupvol/med-bio/aberlang/data.dir/Airwave.dir/chr22_header2_head10000_subset_cols.dose.vcf .')

# Compare samples to test if correct columns were extracted:
cols_extracted <- read.csv('chr22_header2_head10000_subset_cols.dose.vcf', sep = '\t', header = FALSE, 
                           stringsAsFactors = FALSE, strip.white = TRUE, skip = 10)
cols_extracted[1, 1:5]
head(cols_extracted)
dim(cols_extracted)

cols_extracted <- transpose(cols_extracted)
class(cols_extracted)
str(cols_extracted)
head(cols_extracted)
cols_extracted[1:10, 1]
dim(cols_extracted)

merged_headers[1:5, 1]
# Check match with subset of metabolomics samples which have genotyping:
matching <- which(as.character(merged_headers[, 'V1']) %in% as.character(cols_extracted[, 'V1']))
length(matching)

####################


####################
# The end:
# Remove objects that are not necessary to save:

#rm(list=ls(arrayQualityMetrics_raw_cleaned, arrayQualityMetrics_preprocessed))

# To save R workspace with all objects to use at a later time:
# save.image(file=R_session_saved_image, compress='gzip')

#objects_to_save <- (c('normalised_expressed', 'normalised_filtered', 'membership_file_cleaned', 'FAILED_QC_unique'))
#save(list=objects_to_save, file=R_session_saved_image, compress='gzip')
sessionInfo()
q()

# Next: run XGR, IPA, etc, cross with GWAS, ENCODE, etc.
####################
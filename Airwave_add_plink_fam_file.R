###############
# Airwave phenotype data: merge fam file to add sex and other information for plink
# July 2016
# Airwave chronic inflammation proposal
# Objective 

###############


#############################
##Set working directory and file locations and names of required inputs:

# Working directory:
# setwd("/Users/antoniob/Desktop/Airwave_inflammation/results_1.dir")

#Direct output to file as well as printing to screen (plots aren't redirected though, each done separately). 
#Input is not echoed to the output file either.
output_file <- file(paste("R_session_output_Airwave_plink_fam",".txt", sep=""), open='a')
output_file
sink(output_file, append=TRUE, split=TRUE, type = c("output", "message"))

#If the script can run from end to end, use source() with echo to execute and save all input 
#to the output file (and not truncate 150+ character lines):
#source(script_file, echo=TRUE, max.deparse.length=10000)

#Record start of session and locations:
Sys.time()
print(paste('Working directory :', getwd()))

##TO DO extract parameters. This file is very project specific. Check ways of making count comparisons.

# Re-load a previous R session, data and objects:
# load('R_session_saved_image_Airwave_exploratory.RData', verbose=T)

# Filename to save current R session, data and objects at the end:
R_session_saved_image <- paste('R_session_saved_image_plink_fam.RData')

#############################


#############################
## Update packages if necessary and load them:
# source("https://bioconductor.org/biocLite.R")
# biocLite()
#install.packages('foreign')

#############################


#############################
## Load data sets:
covar_data <- read.csv('minimal_covariates_airwave.tsv', header = TRUE, stringsAsFactors = FALSE,
                       strip.white = TRUE, sep = '\t', na.string = c('-9'))

head(covar_data)
file.exists('airwave_b37.fam')
fam_file <-  read.csv('airwave_b37.fam', header = FALSE, stringsAsFactors = FALSE,
                      strip.white = TRUE, sep = ' ')
head(fam_file)

#############################

#############################
# Update fam file with real values from covariates data:
# https://www.cog-genomics.org/plink2/formats#fam
# A text file with no header line, and one line per sample with the following six fields:
# 1. Family ID ('FID')
# 2. Within-family ID ('IID'; cannot be '0')
# 3. Within-family ID of father ('0' if father isn't in dataset)
# 4. Within-family ID of mother ('0' if mother isn't in dataset)
# 5. Sex code ('1' = male, '2' = female, '0' = unknown)
# 6. Phenotype value ('1' = control, '2' = case, '-9'/'0'/non-numeric = missing data if case/control)

head(covar_data)
names(covar_data)
names(fam_file)[1] <- 'FID'
names(fam_file)[2] <- 'IID'

fam_file_2 <- merge(covar_data, fam_file, all.y = TRUE)
head(fam_file_2)
names(fam_file_2)

fam_file_2 <- fam_file_2[, c(1, 2, 12, 13, 11, 10)]
names(fam_file_2)
head(fam_file_2)
tail(fam_file_2)
head(fam_file)
dim(fam_file)
dim(fam_file_2)
dim(covar_data)

#############################


#############################
# Write to file:
write.table(fam_file_2, 'plink_fam_file.tsv', row.names = FALSE, quote = FALSE, sep = '\t', 
            na = '-9', col.names = TRUE)

#############################

#############################
#The end:

# To save R workspace with all objects to use at a later time:
save.image(file=R_session_saved_image, compress='gzip')

q()

# Next: run script for higher level analyses.
#############################
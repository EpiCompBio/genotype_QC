
#################################
#Manhattan plot highlight SNPs of interest

# Requires:
# Input file with data to read with p-values, at the moment only reads Plink's output.
# List of SNPs of interest to highlight in plot
# See:
# http://www.gettinggeneticsdone.com/2011/04/annotated-manhattan-plots-and-qq-plots.html

#Antonio Berlanga-Taylor
#26/01/2016

#################################

args <- commandArgs(trailingOnly = TRUE)
library(qqman)


#Set input and output files:
input_file = as.character(args[1]) #'P140343-Results_FinalReport_clean-individuals_and_SNPs'
SNPs_of_interest = as.character(args[2])

output_path = file.path(getwd(), paste('manhattan_plot_SNPs_of_interest_', 
                                       input_file, '.png', sep = ''))

# Read data:
data = na.omit(read.table(input_file, header = TRUE, na.strings = 'NA'))
#summary(data)
#head(subset(data, select=c(SNP, CHR, BP, P)))

SNPs_of_interest = scan(SNPs_of_interest, character())

#Plot:
png(output_path)
manhattan(data, ylim = c(0, 10), col=c("black","#666666","#CC6600"),
          highlight = SNPs_of_interest, main="Manhattan Plot")
dev.off()

q()


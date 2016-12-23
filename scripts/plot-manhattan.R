#################################
#Manhattan plot

#Requires:
#Input file with data to read with p-values, at the moment only reads Plink's output.


#Antonio Berlanga-Taylor
#26/01/2016

#################################

args <- commandArgs(trailingOnly = TRUE)
library(qqman)


#Set input and output files:
input_file = as.character(args[1]) #'P140343-Results_FinalReport_clean-individuals_and_SNPs'

output_path = file.path(getwd(), paste('manhattan_plot_', input_file, '.png', sep = ''))

# Read data:
data = na.omit(read.table(input_file, header = TRUE, na.strings = 'NA'))
#class(data)
#summary(data)
#head(subset(data, select=c(SNP, CHR, BP, P)))

#Plot:
png(output_path)
manhattan(data, ylim = c(0, 10), col=c("black","#666666","#CC6600"),
          main = 'Manhattan plot')
dev.off()

q()

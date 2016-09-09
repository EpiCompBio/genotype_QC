#################################
#Identity by descent plot
#Plots IBD results from genotype data

#Script modified from www.nature.com/nprot/journal/v5/n9/pdf/nprot.2010.116.pdf

#Requires:
  #Input file with data to read in (.genome file) after plink analysis
  #Set cut-off of IBD, 'IBD_cutoff', currently at 0.1875, typical value used though
  #Runs in current directory

#Antonio Berlanga-Taylor
#05/03/2015

#################################

args1 <- commandArgs(trailingOnly = TRUE)

#Set input and output files:
input_file = as.character(args1) #'P140343-Results_FinalReport'
input_path = file.path(getwd(), input_file)
output_path = file.path(getwd(), 'IBD-plot.pdf')
IBD_cutoff = as.numeric(0.1875)
FAILED_IDs_file = file.path(getwd(), 'fail-IBD-check.FAILED_QC')

#Read in data:
data=read.table(paste(input_path, ".genome", sep=''), h=T)

#Plot:
pdf(output_path)

hist(data$PI_HAT, ylim=c(0,100), col="RED", breaks=100, 
     xlab="Estimated mean pairwise IBD",main="")

abline(v=IBD_cutoff, col="gray32", lty=2)

dev.off()

#Write to table individuals failing IBD results:
out = which(data$PI_HAT > IBD_cutoff)
write.table(data[out,], file = FAILED_IDs_file, append = FALSE, col.names=TRUE, row.names=FALSE, 
            sep='\t', quote= FALSE)

quit()

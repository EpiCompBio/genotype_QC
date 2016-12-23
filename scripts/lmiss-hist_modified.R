#################################
#Variant level missing rates plot
#Plots proportion of missing genotypes variant level and identifies outliers to exclude.

#Script www.nature.com/nprot/journal/v5/n9/pdf/nprot.2010.116.pdf

#Requires:
  #lmiss file from plink output after running: 
  #"plink2 --bfile xxx --missing --out xxx

#Antonio Berlanga-Taylor
#15/01/2016

#################################

args1 <- commandArgs(trailingOnly = TRUE)

#Set input and output files:
input_file = as.character(args1)  #'P140343-Results_FinalReport'

input_path = file.path(getwd(), input_file)

output_path = file.path(getwd(), 'missing_variants_lmiss.png')

# Plot missing variants:
x = read.table(paste(input_path, '.lmiss', sep=''), header = TRUE)

ylabels=c("0", "20K", "40K", "60K", "80K", "100K")
xlabels=c("0.0001", "0.001", "0.01", "0.1", "1")
#par(mfrow=c(1,1))
png(output_path)
hist(log10(x$F_MISS),axes=F,xlim=c(-4,0),col="RED",ylab="Number of SNPs",xlab="Fraction of missing data",main="All SNPs",ylim=c(0,100000))
axis(side=2,labels=F)
mtext(ylabels,side=2,las=2, at=c(0,20000,40000,60000,80000,100000),line=1)
axis(side=1,labels=F)
mtext(xlabels,side=1,at=c(-4,-3,-2,-1,0),line=1)
abline(v=log10(0.05),lty=2)
dev.off()
q()


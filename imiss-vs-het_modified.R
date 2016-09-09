#################################
#Missing data rates and outlying heterozygosity rate plot
#Plots Heterozygosity rate and proportion of missing genotypes, and identifies outliers to exclude.

#Script www.nature.com/nprot/journal/v5/n9/pdf/nprot.2010.116.pdf

#Requires:
  #heterozygosity rate SD (currently at 3)
  #Threshold for proportion of missing genotypes (currently at 0.03)
  #File with data to read in (e.g. Plink's '.imiss' file after running 
  #"plink --bfile raw-GWA-data --missing --out raw-GWA-data"
  #Runs in current directory

#Antonio Berlanga-Taylor
#04/03/2015

#################################

library("geneplotter")

args1 <- commandArgs(trailingOnly = TRUE)

#Set input and output files:
input_file = as.character(args1)  #'P140343-Results_FinalReport'

input_path = file.path(getwd(), input_file)

output_path = file.path(getwd(), 'missing_genotypes_vs_het_rate.pdf')

#Set ablines standard deviation for heterozygosity rate and for proportion of missing genotypes:
het_rate_line_sd = as.numeric(3)

missing_genotypes_proportion = as.numeric(0.03)

print(paste("standard deviation for heterozygosity rate =", het_rate_line_sd, 
            "and proportion of missing genotypes =", missing_genotypes_proportion))

#Read data:
imiss = read.table(paste(input_path, ".imiss", sep=""), h=T)

imiss$logF_MISS = log10(imiss[,6])

het = read.table(paste(input_path, ".het", sep=""), h=T)

het$meanHet = (het$N.NM. - het$O.HOM.) / het$N.NM.

#Plot:
colors <- densCols(imiss$logF_MISS, het$meanHet)

pdf(output_path)

plot(imiss$logF_MISS, het$meanHet, col=colors, xlim=c(-3,0),
     ylim=c(0,0.5), pch=20, xlab="Proportion of missing genotypes",
     ylab="Heterozygosity rate", axes=F)

axis(2, at=c(0,0.05,0.10,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5), tick=T)

axis(1, at=c(-3,-2,-1,0), labels=c(0.001,0.01,0.1,1))

abline(h=mean(het$meanHet) - (het_rate_line_sd*sd(het$meanHet)), col="RED", lty=2)

abline(h=mean(het$meanHet) + (het_rate_line_sd*sd(het$meanHet)), col="RED", lty=2)

abline(v=log10(missing_genotypes_proportion), col="RED", lty=2)

dev.off()

#Print which individuals are outside of het_rate_line_sd and missing_genotypes_proportion:
print(paste('Mean heterozygosity rate is ', round(mean(het$meanHet), 4)))

print(paste('Standard deviation of heterozygosity rate is ', round(sd(het$meanHet), 4)))

het_rate_outliers_above <- het$meanHet > (mean(het$meanHet) + (het_rate_line_sd*sd(het$meanHet)))

which_het_rate_outliers_above <- which(het_rate_outliers_above == TRUE)

print(paste("FAILED QC: FID and IID of individuals above standard deviation for heterozygosity rate 
            (SD>", het_rate_line_sd, "):", het$FID[which_het_rate_outliers_above], 
            het$IID[which_het_rate_outliers_above]))

#head(het)
#het$meanHet[c(64,197)]
#head(het[,c('FID','meanHet')])
#het$meanHet[which_het_rate_outliers_above]

het_rate_outliers_below <- het$meanHet < (mean(het$meanHet) - (het_rate_line_sd*sd(het$meanHet)))

which_het_rate_outliers_below <- which(het_rate_outliers_below == TRUE)

print(paste("FAILED QC: FID and IID of individuals below standard deviation for heterozygosity rate (SD<", 
            het_rate_line_sd, "):", het$FID[which_het_rate_outliers_below], 
            het$IID[which_het_rate_outliers_below]))

missing_genotypes_outliers = imiss$logF_MISS > (log10(missing_genotypes_proportion))

which_missing_genotypes_outliers <- which(missing_genotypes_outliers == TRUE)

print(paste("FAILED QC: FID and IID of individuals above threshold for proportion of missing genotypes (>", 
            missing_genotypes_proportion, "):", imiss$FID[which_missing_genotypes_outliers], 
            imiss$IID[which_missing_genotypes_outliers]))

#head(imiss)
#which(imiss$F_MISS > missing_genotypes_proportion)


#Output failed QC samples to tab delimited table:
#rep('Above standard deviation for heterozygosity rate', byrow=TRUE):
t1 <- matrix(c(het$FID[which_het_rate_outliers_above], het$IID[which_het_rate_outliers_above], 
               round(het$meanHet[which_het_rate_outliers_above], 4)), ncol=3)

#'Below standard deviation for heterozygosity rate':
t2 <- matrix(c(het$FID[which_het_rate_outliers_below], het$IID[which_het_rate_outliers_below],
               round(het$meanHet[which_het_rate_outliers_below], 4)), ncol=3)

#'Above threshold for proportion of missing genotypes':
t3 <- matrix(c(imiss$FID[which_missing_genotypes_outliers], imiss$IID[which_missing_genotypes_outliers],
               round(imiss$F_MISS[which_missing_genotypes_outliers], 4)), ncol=3)

#Convert matrices to table and save to file:
table <- data.frame(rbind(t1, t2, t3))

colnames(table) <- c('FID', 'IID', 'Value') 

write.table(table, file = file.path(getwd(),'missing_genotypes_and_het_rate.FAILED_QC'),
            append = FALSE, col.names = TRUE, quote = FALSE, row.names = FALSE, sep = '\t')

q()

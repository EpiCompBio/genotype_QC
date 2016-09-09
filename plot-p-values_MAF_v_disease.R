#################################
#P-value distribution plot

#Plots distribution of p-values from SNP association analysis after plink --assoc output

#Requires:
#Input file with data to read in after plink analysis
#Set abline for cut-off of p-value, 'p_cutoff', currently at 0.05
#Runs in current directory

#Antonio Berlanga-Taylor
#06/03/2015

#################################

args <- commandArgs(trailingOnly = TRUE)


#Set input and output files:
input_file = as.character(args[1]) #'P140343-Results_FinalReport_clean-individuals_and_SNPs'
#input_path = file.path(getwd(), input_file)
test_done = as.character(args[2])

output_path = file.path(getwd(), paste('p-values-distribution-plot-MAF_v_disease_', test_done, '.png', sep = ''))
p_cutoff = as.numeric(0.05)
output_table = file.path(getwd(), paste('SNPs_under_cut-off_no_adj_MAF_v_disease_', test_done, '.txt', sep = ''))

#Read data:
#data = read.table(paste(input_path, ".assoc", sep=''), h=T)
data = read.table(input_file, header = TRUE)

#head(data)
#dim(data)

#Plot:
png(output_path)
hist(data$P, axes=T, col="BLUE", ylab="Number of SNPs",
     xlab="p-value", main="Distribution of p-values")
abline(v=p_cutoff, col='RED', lty='solid')

dev.off()

#Identify SNPs which are below the p-value cut-off:
out = which(data$P < p_cutoff)

#Write SNPs to table:
write.table(data[out,], file = output_table, append = FALSE, col.names=TRUE,
            row.names=FALSE, sep='\t', quote= FALSE)

#head(data[out,])

q()

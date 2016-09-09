#################################
#Search for SNP IDs in plink output

#Searches for specific SNP IDs after plink analyses in .model output

#Requires:
#Input file .model results with association test analysis
#List of SNPs to search
#P-value cut-off (currently set at 0.05, unadjusted, as output by plink)
#
#Runs in current directory
#Outputs full list and list under cut-off specified

#Antonio Berlanga-Taylor
#06/03/2015

#################################

#Set command line arguments:
args <- commandArgs(trailingOnly = TRUE)

#Set input and output files:
input_file = args[1]
#input_file = 'P140343-Results_FinalReport_clean-individuals_and_SNPs_vitd0.qassoc'

input_path = file.path(getwd(), input_file)
p_cutoff = as.numeric(0.05)

output_file = (paste(getwd(), '/SNP_subset_test_of_assoc_', args[1], sep=''))
#output_file = (paste(getwd(), '/SNP_subset_test_of_assoc_', input_file, sep=''))


output_file_cutoff = (paste(getwd(), '/SNP_subset_test_of_assoc_under_cutoff_unadj_', 
                            args[1], sep=''))
#output_file_cutoff = (paste(getwd(), '/SNP_subset_test_of_assoc_under_cutoff_unadj_', 
#                            input_file, sep=''))


#Read data:
#data = read.table(paste(input_path, ".model", sep=''), h=T)
data = read.table(input_path, h=T)
#input_SNPs = c('rs2282679', 'rs12785878', 'rs10741657', 'rs2282679', 
#  'rs3829251', 'rs11234027', 'rs2060793', 'rs1993116')
input_SNPs = read.table('VD_SNPs_GWAS_list.txt', header = FALSE, sep = '\n')

#head(data)
#dim(data)

#Search for elements:
SNPs_of_interest = input_SNPs

#Intersect with original table:
out = sort(match(SNPs_of_interest$V1, data$SNP))
#out = match(SNPs_of_interest$V1, data$SNP)
out_table = data.frame(data[out, ])

#Subset SNPs under p-value cutoff:
out_p_cutoff = which(out_table$P < p_cutoff)
out_p_cutoff_table = data.frame(out_table[out_p_cutoff, ])

#Write SNPs to files:
write.table(out_table, file = output_file, append = FALSE, col.names=TRUE,
            row.names=FALSE, sep='\t', quote= FALSE)

write.table(out_p_cutoff_table, file = output_file_cutoff, append = FALSE, col.names=TRUE,
            row.names=FALSE, sep='\t', quote= FALSE)

quit()

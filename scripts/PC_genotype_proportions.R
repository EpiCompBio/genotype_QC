setwd('/ifs/projects/proj043/analysis.dir/association_analysis_1.dir')
genotype_PCs <- read.csv('P140343-Results_FinalReport_clean-individuals_and_SNPs_pca.eigenval', sep = ' ', 
                         check.names = FALSE, header = FALSE)
genotype_PCs
# genotype_PCs <- genotype_PCs[, -1]
# rownames(genotype_PCs) <- genotype_PCs[, 1]
# genotype_PCs[, 1] <- NULL
# genotype_PCs[1:5, 1:5]
summary(genotype_PCs)

genotype_PCs_prop$total <- colSums(genotype_PCs)
genotype_PCs_prop$variance <- rowSums(genotype_PCs)/colSums(genotype_PCs)
genotype_PCs_prop$total_var <- sum(genotype_PCs_prop$variance)
genotype_PCs_prop
qplot(genotype_PCs_prop$variance)

genotype_PCs <- data.frame(table(genotype_PCs[]))
genotype_PCs$proportion <- prop.table(genotype_PCs$Frequency)
genotype_PCs$cum_prop <- cumsum(genotype_PCs$proportion)


mkMyTable <- function(X){
  Table <- data.frame( table(X) )
  Table$Prop <- prop.table( Table$Freq )
  Table$CumProp <-  cumsum( Table$Prop )
  Table
}

myTable <- mkMyTable(genotype_PCs[, 3:12]) 
#########################

genotype_data <- read.csv('P140343-Results_FinalReport_clean-individuals_and_SNPs', sep = ' ', 
                         check.names = FALSE, header = FALSE)
pca_normalised_filtered <- prcomp(t(normalised_filtered$E), center=TRUE, scale=TRUE)

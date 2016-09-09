#################################
#Ancestry divergence PCA plot
#Plots a scatter diagram using the first principal components after running smartpca.perl

#Script modified from www.nature.com/nprot/journal/v5/n9/pdf/nprot.2010.116.pdf

#Requires:
  #Input file with data to read in after plink and smartpca.perl analysis
  #Set abline for cut-off of divergent ancestry, 'divergence_cutoff', currently at 
  #0.072 of PC2 (arbitrary)
  #Runs in current directory

#Antonio Berlanga-Taylor
#05/03/2015

#################################
args1 <- commandArgs(trailingOnly = TRUE)

#Set input and output files:
input_file = as.character(args1)  #'P140343-Results_FinalReport.Anderson_hapmap_snps_'
input_path = file.path(getwd(), input_file)
output_path = file.path(getwd(), 'pca-ancestry-plot.pdf')
divergence_cutoff = as.numeric(0.072)
FAILED_IDs_file = file.path(getwd(), 'fail-ancestry_divergence-check.FAILED_QC')

#Read data:
data=read.table(paste(input_path, "hapmap3r2.pruned.pca.evec", sep=''), h=F, skip=1, 
                colClasses=c('character', 'numeric', 'numeric', 'character'), 
                col.names=c('ID', 'PC1', 'PC2', 'Group'))

head(data$Group)

#Set variables:
cont=which(data$Group=="Control")
case=which(data$Group=="Case")
CEU=which(data$Group=="3")
CHB=which(data$Group=="4")
JPT=which(data$Group=="5")
YRI=which(data$Group=="6")

#Plot PCA results:
pdf(output_path)
plot(0,0,pch="",xlim=c(-0.1,0.05),ylim=c(-0.05,0.1),xlab="Principal component 1", 
     ylab="Principal component 2")

#Set colours for populations:
points(data$PC1[JPT],data$PC2[JPT],pch=20,col="PURPLE")
points(data$PC1[CHB],data$PC2[CHB],pch=20,col="YELLOW")
points(data$PC1[YRI],data$PC2[YRI],pch=20,col="GREEN")
points(data$PC1[CEU],data$PC2[CEU],pch=20,col="RED")

#Controls are '*' and blue, cases are '+' and black :
points(data$PC1[cont],data$PC2[cont],pch="*",col="BLUE")
points(data$PC1[case],data$PC2[case],pch="+",col="BLACK")
par(cex=2)

#Set cut-off line, currently at 0.072 of PC2 (arbitrary):
abline(h=divergence_cutoff, col="gray32", lty=2)

dev.off()

#Identify individuals which are outliers and are to be excluded:
out = which(data$PC2 < divergence_cutoff)
Cases = which((data$Group == 'Case') & data$PC2 < divergence_cutoff)
Controls = which((data$Group == 'Control') & data$PC2 < divergence_cutoff)

out_study = c(Cases, Controls)

#Separate ID column from input to create FID and IID columns 
#(to grep out later for exclusion by plink):

data$FID = as.character(lapply(strsplit(as.character(data[[1]]), split=":"), '[', 1))
data$IID = as.character(sapply(strsplit(as.character(data[[1]]), split=":"), '[', 2))

#Reorder columns to have FID and IID written first:

data_reordered = data[out_study,][c(5,6,1:4)]

#Write to table individuals from study which fail the ancestry divergence QC:
write.table(data_reordered, file = FAILED_IDs_file, append = FALSE, 
            row.names=FALSE, sep='\t', quote= FALSE)



#head(data)
#tail(data$ID)
#dim(data)
#sum(data$PC2)
#max(data$PC2)

quit()

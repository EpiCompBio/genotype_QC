#!/bin/bash

#########################
#Genotype QC pipeline
#Removal of low-quality SNP markers
#########################
#Antonio J Berlanga-Taylor
#14 January 2016

#########################
##Quality control of subject and SNP genotyping data
# Usually run individuals QC first and then markers QC (this file). If prioritising keeping individuals over markers for a study with low sample size then
# run marker QC first. 
# Note this is now using plink 1.9, which is currently in beta

#Check:
#Anderson et al. 2010 protocol: http://www.nature.com/nprot/journal/v5/n9/pdf/nprot.2010.116.pdf
#Winkler et al. 2014 protocol (meta-analysis of GWAS): http://www.nature.com/nprot/journal/v9/n5/pdf/nprot.2014.071.pdf
#Plink tutorial: http://pngu.mgh.harvard.edu/~purcell/plink/tutorial.shtml

# Inputs are bed, bim, and fam, and outputs are plink files with clean marker set which passed QC.

#Dependencies:
#EIGENSOFT: ftp://pricelab:pricelab@ftp.broadinstitute.org/EIGENSOFT/EIG6.0.1.tar.gz
#Plink: http://pngu.mgh.harvard.edu/~purcell/plink/ibdibs.shtml
#Plotting scripts and others from Anderson et al 2010 protocol, currently in:
#/ifs/projects/proj043/analysis.dir/genotypes.dir
#and
#/ifs/projects/proj043/analysis.dir/genotypes.dir/Anderson_et_al_2010_Nat_Protocols.dir/

##Plink help:
#plink2 --noweb --help

## TO DO:
# Modify to run in ruffus...

#########################


#########################
# Set files up:
ln -s /ifs/projects/proj043/analysis.dir/genotypes.dir/P140343-Results_FinalReport.bim .
ln -s /ifs/projects/proj043/analysis.dir/genotypes.dir/P140343-Results_FinalReport.bed .
ln -s /ifs/projects/proj043/analysis.dir/genotypes.dir/P140343-Results_FinalReport.fam .

# Generate a baseline report of allele frequencies:
plink2 --bfile P140343-Results_FinalReport --freqx --out P140343-Results_FinalReport


#1)	Identification of SNPs with excessive missing call rates
#Calculate the missing genotype rate for each marker. Output will be .lmiss (variant-based missing data report) and .imiss (sample-based missing data report):

plink2 --noweb --bfile P140343-Results_FinalReport --missing --out P140343-Results_FinalReport

#Plot a histogram of the missing genotype rate (.lmiss data) to identify a threshold for extreme genotype failure rate. Use column 5 of the .lmiss file 
#and the lmiss-hist.Rscript file. A call-rate threshold of 3% is suggested.
# Run R script with file (minus suffix as argument)

Rscript lmiss-hist_modified.R P140343-Results_FinalReport

#2)	Identification of differing genotype call rates between cases and controls
#Test all markers for differences in call rate between cases and controls. Output will be in .missing file:

plink2 --noweb --bfile P140343-Results_FinalReport --test-missing --out P140343-Results_FinalReport

#Create a file that contains all SNPs with a significantly different (P<0.00001) missing data rate between cases and controls. 
#Output will be in 'fail-diffmiss-qc.txt'.
# I modified the file run-diffmiss-qc.pl to run-diffmiss-qc_p-value_0.01.pl with higher p-value:
 
perl run-diffmiss-qc_p-value_0.01.pl P140343-Results_FinalReport

#3)	SNP quality (filtering of monomorphic SNPs, SNPs with missing values or nonsense values; imputation quality; low call rate; 
#violation of Hardy-Weinberg Equilibrium; duplication; and minimum allele frequency). 
# Exclude all non-autosomal and unplaced variants including mitochondrial, X and Y chromosomes (with --autosome or use --chr and specify, eg --chr 1-22):
#Remove poor SNPs and create a clean GWA data file. SNPs failing previous QC steps, those with an MAF<0.10 and an HWE P-value < 0.00001 (in controls) 
#are removed. 

#Filtering criteria of SNPs for BEST-D:
# SNPs with significantly different missing rate between cases and controls (fail-diffmiss file from steps 1 to 3)
# Sample and SNP call rates above 98%
# Minor allele frequency (MAF) of 10% 
# Hardy-Weinberg equilibrium threshold of 1x10-6
# Markers which are unplaced or not in chr 1-22 (all non-autosomal SNPs)

plink2 --noweb --bfile P140343-Results_FinalReport --maf 0.10 --geno 0.02 --hwe 0.00001 --autosome --exclude fail-diffmiss-qc.txt --make-bed \
	--out P140343-Results_FinalReport_clean_SNPs_autosome

# Generate a report of allele frequencies after removal of low quality markers:
plink2 --bfile P140343-Results_FinalReport_clean_SNPs_autosome --freqx --out P140343-Results_FinalReport_clean_SNPs_autosome



#Individuals and SNPs failing QC steps will be removed and reported at each step. 
#References: PLINK (Purcell, Neale et al., 2007), related software and protocols (Anderson, Pettersson et al., 2010; Winkler, Day et al., 2014; Sham and #Purcell, 2014).

#TO DO / CHECK if needed:
#High LD regions removal (get file)
#Look for finer population substructure

#########################


#########################
#Next steps:
#Individuals QC or Association test analysis, see separate files.

#########################

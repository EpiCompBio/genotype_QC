#!/bin/bash

# Run without set -e as plink errors with unflipped SNPs, when it does it saves the problem SNPs to '.missnp' file and the next command can then run using that file.

#set -e 

#########################
#Genotype QC pipeline
#Removal of low-quality individuals from genotype file
#########################
#Antonio J Berlanga-Taylor
#14 January 2016

#########################
##Quality control of subject and SNP genotyping data
# Note this is now using plink 1.9, which is currently in beta

#Check:
#Anderson et al. 2010 protocol: http://www.nature.com/nprot/journal/v5/n9/pdf/nprot.2010.116.pdf
#Winkler et al. 2014 protocol (meta-analysis of GWAS): http://www.nature.com/nprot/journal/v9/n5/pdf/nprot.2014.071.pdf
#Plink tutorial: http://pngu.mgh.harvard.edu/~purcell/plink/tutorial.shtml

# Inputs are bed, bim, fam files from plink (after pre-processing and optionally individual marker QC)
# Outputs are clean genotype file for downstream association, eQTL, etc. analysis.

#Dependencies:
#EIGENSOFT: ftp://pricelab:pricelab@ftp.broadinstitute.org/EIGENSOFT/EIG6.0.1.tar.gz
#Plink: http://pngu.mgh.harvard.edu/~purcell/plink/ibdibs.shtml
#Plotting scripts and others from Anderson et al 2010 protocol, currently in:
#/ifs/projects/proj043/analysis.dir/genotypes.dir
#and
#/ifs/projects/proj043/analysis.dir/genotypes.dir/Anderson_et_al_2010_Nat_Protocols.dir/
# Should all be already in /ifs/devel/antoniob/projects/BEST-D/

#########################
bfile=$1
#bfile="P140343-Results_FinalReport"
SNPs_of_interest=$2
#SNPs_of_interest="gwas-association-downloaded_2016-07-28-fibrinogen.tsv"
base_directory_scripts=$3
#base_directory_scripts="/home/aberlang/bin/bin_symlinks/"
#########################
##a.	Gender mis-identification
#See http://www.nature.com/nprot/journal/v5/n9/pdf/nprot.2010.116.pdf (reference as above)

#Calculate the mean homozygosity rate across X-chromosome markers for each individual in the study:
# If running after markers QC (to prioritise keeping individuals over markers for a study with low sample size), run sexcheck on raw genotype file which
# contains non-autosome chromosomes.

plink --bfile $bfile --check-sex --out $bfile

#Produce a list of individuals with discordant sex data:

grep PROBLEM ${bfile}.sexcheck > gender_mis-identification_check.FAILED_QC

#The file will contain the family IDs (column 1) and individual ID (column 2) for these individuals. #Column 3 denotes ascertained sex and #column 4 denotes
# sex according to genotype data (1=male, 2=female). 
# When the homozygosity rate is more than 0.2 but less than 0.8, the genotype data are inconclusive regarding the sex of an individual and these are 
# marked in column 4 with a 0.

#Report the IDs of individuals with discordant sex information to those who conducted sex phenotyping. If discrepancy cannot be resolved, 
# add the family ID (FID) and individual ID (IID) of the samples to a QC file and exclude from further analysis. 

#b.	Individuals with elevated missing data rates or outlying heterozygosity rate; 

#Create the plink files xxxx.imiss and xxxx.lmiss for sample based and variant based missing SNPs, respectively. xxxx.lmiss QC is used in the marker QC.
# Here, the fourth column in the .imiss file (N_MISS) denotes the number of missing SNPs 
#and the sixth column (F_MISS) denotes the proportion of missing SNPs per individual:

plink --bfile ${bfile}_clean_SNPs_autosome --missing --out ${bfile}_clean_SNPs_autosome

#Create the plink file xxxx.het, in which the third column denotes the observed number of homozygous genotypes [O(Hom)] and the fifth column denotes 
#the number of nonmissing genotypes [N(NM)] per individual:

plink --bfile ${bfile}_clean_SNPs_autosome --het --out ${bfile}_clean_SNPs_autosome


#Calculate the observed heterozygosity rate per individual using the formula:
#(N(NM) + O(Hom)) / N(NM)

#Create a graph in which the observed heterozygosity rate per individual is plotted on the x axis and the proportion of missing SNPs per individuals is 
#plotted on the y axis:

#See script /ifs/projects/proj043/analysis.dir/genotypes.dir/imiss-vs-het_modified.R
#Modified from /ifs/projects/proj043/analysis.dir/genotypes.dir/Anderson_et_al_2010_Nat_Protocols.dir/imiss-vs-het.R

#At the cmd line run: 

Rscript ${base_directory_scripts}imiss-vs-het_modified.R ${bfile}_clean_SNPs_autosome

#Examine the plot to decide QC thresholds at which to exclude individuals based on elevated missing genotypes or extreme heterozygosity. 
#Ref. above excluded individuals with a genotype failure rate >= 0.03 (vertical dashed line) and/or a heterozygosity rate ± 3 s.d. from the mean
#(horizontal dashed lines).

#Add FID and IID of samples which fail QC:
#The script imiss-vs-het.R generates a table missing_genotypes_and_het_rate.FAILED_QC with the IDs of individuals who failed this QC.

#c1.	Identification of related or duplicated individuals

#Reduce the number of SNPs used to create the identity by descent (IBS) matrix by pruning the data set so that no pair of SNPs (within a given number of 
#base pairs) has an r2 value greater than a given threshold (typically 0.2). This reduces computational complexity.
#Check http://pngu.mgh.harvard.edu/~purcell/plink/summary.shtml#prune
#First create the list of SNPs that are in LD with each other from the raw data file:

plink --bfile ${bfile}_clean_SNPs_autosome --indep-pairwise 50 5 0.2

#The command specifies 50 5 0.2 which would a) consider a window of 50 SNPs, b) calculate LD between each pair of SNPs in the window, 
#c) remove one of a pair of SNPs if the LD is greater than 0.2, d) shift the window 5 SNPs forward and repeat the procedure. 
#--indep-pairwise parameters are window size, step and r^2 threshold. Another way of filtering is with sliding windows instead of pairwise (--indep).

#Two files are generated, plink.prune.in and plink.prune.out. 
#Generate pairwise IBS for all pairs of individuals based on the reduced marker set:

plink --bfile ${bfile}_clean_SNPs_autosome --extract plink.prune.in --genome \
	--out ${bfile}_clean_SNPs_autosome_pruned

#With plink 1.07 I previously ran this with qsub, files are in /ifs/projects/proj043/analysis.dir/genotypes.dir/qsub_params_for_plink_IBS.sh 
#and
#qsub_plink_IBS.sh

#qsub qsub_params_for_plink_IBS.sh

#This generates a .genome file amongst others. 

#Change file names so IBD R script can run (.imiss and .genome must have same prefix):

ln -s ${bfile}_clean_SNPs_autosome_pruned.genome ${bfile}_clean_SNPs_autosome.genome

#To identify all pairs of individuals with an IBD > 0.1875, run R script below (could also use perl run-IBD-QC.pl but R results more complete 
#plus generates a plot):

Rscript ${base_directory_scripts}plot-IBD_modified.R ${bfile}_clean_SNPs_autosome

#IBD cut-offs Anderson 2010 explanation: 
#IBD>0.98 for duplicates or monozygotic twins, IBD=0.5 for first-degree relatives, IBD=0.25 for second-degree relatives and
#IBD=0.125 for third-degree relatives. Owing to genotyping error, LD and population structure, there is often some variation.
#It is typical to remove one individual from each pair with an IBD value of >0.1875, halfway between third- and second-degree relatives. 


#Perl script is from Anderson 2010 Nat Protocols, in /ifs/projects/proj043/analysis.dir/genotypes.dir/Anderson_et_al_2010_Nat_Protocols.dir/run-IBD-QC.pl
#This looks at the individual call rates stored in xxx.imiss and outputs the IDs of the individual with the lowest call rate to 'fail-IBD-QC xxx' 
#for subsequent removal.


#c2.	Identification of individuals of divergent ancestry

#Merge study genotypes to HapMap Phase III (HapMap3) data from four ethnic populations. 
#I'm currently using SNPs from Anderson et al. 2010. 
#TO DO: Finer population sub-structure requires more populations (e.g. within Europe).
# (a) Create a new BED file, excluding SNPs from the study data that do not feature in the genotype data of the four original HapMap3 populations:

#ln -s /ifs/projects/proj043/analysis.dir/genotypes.dir/Anderson_et_al_2010_Nat_Protocols.dir/hapmap3r2_CEU.CHB.JPT.YRI.* .

plink --bfile ${bfile}_clean_SNPs_autosome --extract hapmap3r2_CEU.CHB.JPT.YRI.no-at-cg-snps.txt --make-bed \
	--out ${bfile}_clean_SNPs_autosome.Anderson_hapmap_snps


# (b) Merge the raw-GWA-data.hapmap-snps files with the HapMap data and extract the pruned SNP set:
# This errors due to +/- strand problem, use --flip as below.

plink --bfile ${bfile}_clean_SNPs_autosome.Anderson_hapmap_snps \
	--bmerge hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps \
	--extract plink.prune.in --make-bed --out ${bfile}_clean_SNPs_autosome.Anderson_hapmap_snps_hapmap3r2.pruned

#If this errors with (in plink 1.07) ´ERROR: Stopping due to mis-matching SNPs -- check +/- strand?´, it is likely due to problem
# with strand flipping. 
#The alleles at each marker must be aligned to the same DNA strand to allow the study data to merge correctly. 
#Because not all SNPs are required for this analysis, A->T and C->G SNPs, which are more difficult to align, can be omitted. 
#Problem SNPs are saved by plink to the .missnp file. 
#Re-run command above (a) with --flip .pruned-merge.missnp:

plink --bfile ${bfile}_clean_SNPs_autosome --extract hapmap3r2_CEU.CHB.JPT.YRI.no-at-cg-snps.txt \
	--flip ${bfile}_clean_SNPs_autosome.Anderson_hapmap_snps_hapmap3r2.pruned-merge.missnp \
	--make-bed \
	--out ${bfile}_clean_SNPs_autosome.Anderson_hapmap_snps

#Then re-run (b) (the --bmerge command) exactly as above:
plink --bfile ${bfile}_clean_SNPs_autosome.Anderson_hapmap_snps --bmerge hapmap3r2_CEU.CHB.JPT.YRI.founders.no-at-cg-snps \
	--extract plink.prune.in --make-bed --out ${bfile}_clean_SNPs_autosome.Anderson_hapmap_snps_hapmap3r2.pruned


#Create a copy of the output BIM and FAM files just generated:
cp ${bfile}_clean_SNPs_autosome.Anderson_hapmap_snps_hapmap3r2.pruned.bim \
	${bfile}_clean_SNPs_autosome.Anderson_hapmap_snps_hapmap3r2.pruned.pedsnp

cp ${bfile}_clean_SNPs_autosome.Anderson_hapmap_snps_hapmap3r2.pruned.fam \
	${bfile}_clean_SNPs_autosome.Anderson_hapmap_snps_hapmap3r2.pruned.pedind

#Conduct a PCA on the merged data (requires smartpca.pl from EIGENSOFT, see dependencies above). Files are from the previous output 
#plus 'pca-populations.txt' file that specifies populations columns:
# smartpca is part of ftp://pricelab:pricelab@ftp.broadinstitute.org/EIGENSOFT/EIG-6.1.2.tar.gz
 
#${base_directory_scripts}smartpca.pl -i ${bfile}_clean_SNPs_autosome.Anderson_hapmap_snps_hapmap3r2.pruned.bed \
#        -a ${bfile}_clean_SNPs_autosome.Anderson_hapmap_snps_hapmap3r2.pruned.pedsnp \
#        -b ${bfile}_clean_SNPs_autosome.Anderson_hapmap_snps_hapmap3r2.pruned.pedind \
#        -o ${bfile}_clean_SNPs_autosome.Anderson_hapmap_snps_hapmap3r2.pruned.pca \
#        -p ${bfile}_clean_SNPs_autosome.Anderson_hapmap_snps_hapmap3r2.pruned.plot \
#        -e ${bfile}_clean_SNPs_autosome.Anderson_hapmap_snps_hapmap3r2.pruned.eval \
#        -l ${bfile}_clean_SNPs_autosome.Anderson_hapmap_snps_hapmap3r2.pruned.log \
#        -k 2 -t 2 -w pca-populations.txt

# Updated to use flashpca instead (https://github.com/gabraham/flashpca):
# Running as suggested but repeating some steps:
# TO DO: clean up update to flashpca and remove steps used before for smartpca

# Prune data based on LD and exclude regions:
plink --bfile ${bfile}_clean_SNPs_autosome --indep-pairwise 1000 50 0.05 --exclude range exclusion_regions_hg19.txt --out flashpca

# Exract variants from the pruning and create new file (10,000 to 50,000 SNPs are needed):
plink --bfile ${bfile}_clean_SNPs_autosome --extract flashpca.prune.in --make-bed --out ${bfile}_clean_SNPs_autosome_pruned

# Run flashpca:
${base_directory_scripts}flashpca_x86-64 --mem low --bfile ${bfile}_clean_SNPs_autosome_pruned --numthreads 8 --suffix _${bfile}_clean_SNPs_autosome_pruned.txt

#Create a scatter diagram of the first two principal components, including all individuals in the file raw-GWA-data.hapmap3r2.pruned.pca.evec 
#(the first and second principal components are columns 2 and 3, respectively). Use the data in column 4 to color the points according to sample origin. 
#Use the modified Anderson et al. 2010 plot-pca-results_modified.R to plot and extract IDs from individuals failing this QC:
# This script is for smartpca:

Rscript ${base_directory_scripts}plot-pca-results_modified.R ${bfile}_clean_SNPs_autosome.Anderson_hapmap_snps_

#Check the thresholds to apply for exclusion of individuals based on ancestry. Currently at 0.072 for PC2 as suggested.

# Run for flashpca results:


#TO DO: Run fine-scale ancestry stratification analysis with closer references within Europe using CEU, TSI, GBR, FIN and IBS samples from 
#www.1000genomes.org. Robust identification of fine-scale population structure often requires the construction of many (2-10) 
#principal components (Anderson et al 2010 protocol suggestions).

#	Removal of all individuals failing the QC checks
#Concatenate all files listing individuals who fail the previous QC steps into a single file with unique IDs:

cat *FAILED_QC | sort -k1 | head -n -3 > FAILED-QC-INDIVIDUALS.txt

#Remove individuals who failed QC from the data set:

plink --bfile ${bfile}_clean_SNPs_autosome --remove FAILED-QC-INDIVIDUALS.txt --make-bed \
	--out ${bfile}_clean_SNPs_autosome_individuals

#Individuals and SNPs failing QC steps will be removed and reported at each step. 

# Generate a report of allele frequencies after removal of low quality individuals:

plink --bfile ${bfile}_clean_SNPs_autosome_individuals --freq --out ${bfile}_clean_SNPs_autosome_individuals
plink --bfile ${bfile}_clean_SNPs_autosome_individuals --freqx --out ${bfile}_clean_SNPs_autosome_individuals

# Get frequency report for SNPs of interest:

grep -wf $SNPs_of_interest ${bfile}_clean_SNPs_autosome_individuals.frq > \
        ${SNPs_of_interest}.frq

grep -wf $SNPs_of_interest ${bfile}_clean_SNPs_autosome_individuals.frqx > \
	${SNPs_of_interest}.frqx

#References: PLINK (Purcell, Neale et al., 2007), related software and protocols (Anderson, Pettersson et al., 2010; Winkler, Day et al., 
#2014; Sham and #Purcell, 2014).

#########################

#########################
#Next steps:
#Marker QC if not done and association test analysis, see separate files.

#########################

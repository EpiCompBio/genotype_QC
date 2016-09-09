#!/bin/bash

set -e
set -u 

#########################
#Genotype association analysis pipeline
#########################
#Antonio J Berlanga-Taylor
#05 March 2015
# Updated 18 January 2016

#########################
##Association analysis of SNP genotyping data

#Check:
#Anderson et al. 2010 protocol: http://www.nature.com/nprot/journal/v5/n9/pdf/nprot.2010.116.pdf
#Clarke et al. 2011 protocol: http://www.nature.com/nprot/journal/v6/n2/pdf/nprot.2010.182.pdf
#Plink tutorial: http://pngu.mgh.harvard.edu/~purcell/plink/tutorial.shtml


#Dependencies:
#Plink: http://pngu.mgh.harvard.edu/~purcell/plink/ibdibs.shtml
#Haploview: http://www.broadinstitute.org/haploview/haploview
#SNPSpD, calculates the effective number of independent SNPs among a collection of SNPs in LD with each other: http://genepi.qimr.edu.au/general/daleN/SNPSpD/
#Plotting scripts and others currently in: /ifs/devel/antoniob/projects/cgat_projects/BEST-D/genotype_analysis/

#Requirements (generated with Plink or other tool):
#Standard FAM, PED and MAP or binary files after QC analyses (markers QC and individuals QC, see /ifs/devel/antoniob/projects/cgat_projects/BEST-D/genotype_analysis/ files).
#Alternate phenotype and covariate files if available


#########################
#Genotype association test analysis

# See in Data Analysis Protocol for BEST-D:
# "Association tests will be carried out as described previously (Sham and Purcell, 2014; Pirinen, Donnelly et al., 2013) using 
#frequentist methods (typically based on Chi-square tests) and allowing for multiple testing correction. We will consider Bayesian #methods for specific 
#instances, such as fine mapping, as previously described (Stahl, Wegmann et al., 2012; Stephens and Balding, #2009). 
#Association tests will be carried out at the appropriate step for the particular analysis (e.g. during QTL analyses)."
# etc.

#Overview of steps:

#Descriptive summary
#Candidate SNP association test
#Data visualisation
#Multiple testing adjustment
#Correction for population stratification


#########################
##Set-up files and directory

#Output from QC is in: /ifs/projects/proj043/analysis.dir/genotypes.dir
#Clean file after SNP and individual outlier detection is: P140343-Results_FinalReport_clean_SNPs_autosome_individuals
#Directory for analysis is: /ifs/projects/proj043/analysis.dir/

#Create directory:
#mkdir association_analysis_1.dir
#cd association_analysis_1.dir

# Get input files with softlinks:
ln -s /ifs/projects/proj043/analysis.dir/genotypes_2.dir/P140343-Results_FinalReport_clean_SNPs_autosome_individuals.bed .
ln -s /ifs/projects/proj043/analysis.dir/genotypes_2.dir/P140343-Results_FinalReport_clean_SNPs_autosome_individuals.bim .
ln -s /ifs/projects/proj043/analysis.dir/genotypes_2.dir/P140343-Results_FinalReport_clean_SNPs_autosome_individuals.fam .

# Get phenotype file: 
#ln -s /ifs/projects/proj043/analysis.dir/BEST-D_alternate_phenotypes_kit_ID_twice_recoded_Txsubgroups_new_variables.txt .
ln -s /ifs/projects/proj043/analysis.dir/BEST-D_alternate_phenotypes_kit_ID_twice_recoded_Txsubgroups_new_variables.txt final_phenotype_data.tab

# File to pass with covariates that will need to be accounted for in regression:
ln -s /ifs/projects/proj043/analysis.dir/covariates_list_adjustment.txt .

# File with SNPs of interest:
#ln -s /ifs/projects/proj043/analysis.dir/VD_SNPs_GWAS_list.txt .

# Scripts needed for processing after association analysis:
#ln -s /ifs/devel/antoniob/projects/cgat_projects/BEST-D/genotype_analysis/BEST-D_association_analysis_commands.sh .
#ln -s /ifs/devel/antoniob/projects/cgat_projects/BEST-D/genotype_analysis/qsub_BEST-D_association_analysis_commands.sh .
ln -s /ifs/devel/antoniob/projects/cgat_projects/BEST-D/genotype_analysis/plot-p-values_* .
ln -s /ifs/devel/antoniob/projects/cgat_projects/BEST-D/genotype_analysis/extract_SNPs_assoc_test.R .
ln -s /ifs/devel/antoniob/projects/cgat_projects/BEST-D/genotype_analysis/plink_output_processing.R .
ln -s /ifs/devel/antoniob/projects/cgat_projects/BEST-D/genotype_analysis/plot-qqplot.R .
ln -s /ifs/devel/antoniob/projects/cgat_projects/BEST-D/genotype_analysis/plot-manhattan* .

#########################


#########################
#Descriptive summary
#Obtain a summary of MAFs in case and control populations and an estimate of the OR for association between the minor allele (based on the whole sample) 
#and disease:

/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals --assoc \
	--out P140343-Results_FinalReport_clean_SNPs_autosome_individuals

# or run --assoc instead of --model
#Output file is 'xxx.assoc/model'. It has one row per SNP containing the chromosome [CHR], the SNP identifier [SNP], the base-pair location [BP], 
#the minor allele [A1], the frequency of the minor allele in the cases [F_A] and controls [F_U], the major allele [A2] and statistical data for an 
#allelic association test including the chi2-test statistic [CHISQ], the asymptotic P value [P] and the estimated OR for association between the minor 
#allele and disease [OR].

# Plink manual:
# Given a case/control phenotype, --assoc writes the results of a 1df chi-square allelic test to plink.assoc (or .assoc.fisher with 'fisher'/'fisher-midp'), 
# while --model performs four other tests as well (1df dominant gene action, 1df recessive gene action, 2df genotypic, and Cochran-Armitage trend) 
# and writes the combined results to plink.model.
# See candidate SNP association tests below.

#Sanity plot of p-value distribution and output of unadjusted SNPs under specified cut-off, specify filename and test used (assoc or model):

/ifs/apps/apps/R-3.2.3/bin/Rscript plot-p-values_MAF_v_disease.R P140343-Results_FinalReport_clean_SNPs_autosome_individuals.assoc assoc

#########################


#########################
### Candidate SNP association test
#See http://pngu.mgh.harvard.edu/~purcell/plink/anal.shtml for the different analysis options and meaning: 
#--assoc, --model, --fisher, --linear and --logistic 

## Without covariates:

/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals --model \
	--out P140343-Results_FinalReport_clean_SNPs_autosome_individuals

# This also runs the basic descriptive statistics as above (with --assoc or --model). Case/control definitions come from FAM file, change these accordingly.

#The --cell option allows the minimum number of observations in each cell for chi-square tests (default is 5). 
#Output in in the 'data.model' file containing chromosome [CHR], SNP identifier [SNP], minor allele [A1], major allele [A2], 
#test performed [TEST: GENO (genotypic association); TREND (Cochran-Armitage trend); ALLELIC (allelic association); 
#DOM (dominant model); and REC (recessive model)], the cell frequency counts for cases [AFF] and controls [UNAFF], 
#the chi2 test statistic [CHISQ], the degrees of freedom for the test [DF] and the asymptotic P value [P].

#Plot of p-value distribution and output of unadjusted SNPs under specified cut-off:

/ifs/apps/apps/R-3.2.3/bin/Rscript plot-p-values_MAF_v_disease.R P140343-Results_FinalReport_clean_SNPs_autosome_individuals.model model


#Create a file containing output of association tests based on logistic regression assuming a multiplicative model 
#and including covariates in the GWA data. Phenotypes (case/control labels) are in the FAM file. The alternate phenotype file would include covariates
# that can be used for additional adjusting (as well as gender for example).
# Not for BEST-D though

#/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals --logistic sex \
#	--covar BEST-D_alternate_phenotypes_kit_ID_twice_recoded.txt --covar-name bmi0,vitd0 \
#	--out P140343-Results_FinalReport_clean_SNPs_autosome_individuals



#Output is in 'data.assoc.logistic'.
#If no model option is specified, the first row for each SNP corresponds to results for a multiplicative test of association. 
#For the '--genotypic' option, the first row will correspond to a test for additivity and the subsequent row to a separate test for 
#deviation from additivity. 
#For the '--dominant' or '--recessive' models, the first row will correspond to tests for a dominant or recessive model of association, respectively. 
#If covariates were included, each of these P values is adjusted for the effect of the covariates. 
#The subsequent rows for each SNP correspond to separate tests of significance for each of the C covariates included in the regression model. 

#For the '--genotypic' model, there is a final row per SNP corresponding to a 2 d.f. 
#LR test of whether both the additive and the deviation from additivity components of the regression model are significant. 
#Each row contains the chromosome [CHR], the SNP identifier [SNP], the base-pair location [BP], the minor allele [A1], 
#the test performed [TEST: ADD (multiplicative model or genotypic model testing additivity), GENO_2DF (genotypic model), 
#DOMDEV (genotypic model testing deviation from additivity), DOM (dominant model) or REC (recessive model)], the number of missing individuals 
#included [NMISS], the OR, the coefficient z-statistic [STAT] and the asymptotic P value [P].

#Note: ORs for main effects cannot be interpreted directly when interactions are included in the model; 
#their interpretation depends on the exact combination of variables included in the model.


#########################


#########################
#Association test using the alternate phenotype file

#For quantitative traits --assoc uses the Wald test, which will give the same results as --linear by itself. 
# TO DO: for eGFR: check:
# http://www.renal.org/information-resources/the-uk-eckd-guide/about-egfr#sthash.XIjPTFWQ.dpbs
# http://egfrcalc.renal.org/ 

#TO DO: request permutations
# BMI_DEXA12 has many missing variables (32); physical_activity12 (14); so excluding from regression (in 
# 'BEST-D_alternate_phenotypes_kit_ID_twice_recoded_Txsubgroups_new_variables.txt').
#Test one phenotype at a time

#vitd0:
/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals \
	--pheno final_phenotype_data.tab --pheno-name vitd0 \
	--linear sex hide-covar --adjust --ci 0.95 --missing-phenotype -99999 \
	--covar final_phenotype_data.tab \
	--covar-name incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, bmi0, calendar_age_ra, season_randomisation_2 \
	--out QCd_genotypes_vitd0_all_covar_linear \
	--perm 

# vitd0 with --assoc and qt-means but only sex as covariate:
/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals \
	--pheno final_phenotype_data.tab --pheno-name vitd0 \
	--assoc qt-means --adjust --ci 0.95 --missing-phenotype -99999 \
	--out QCd_genotypes_vitd0_all_covar_assoc \
        --perm

#vitd0 with genotypic to test for other than additive effects::
/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals \
	--pheno final_phenotype_data.tab --pheno-name vitd0 \
	--linear sex genotypic hide-covar --adjust --ci 0.95 --missing-phenotype -99999 \
	--covar final_phenotype_data.tab \
	--covar-name incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, bmi0, calendar_age_ra, season_randomisation_2 \
	--out QCd_genotypes_vitd0_genotypic_all_covar \
        --perm

#To specify a genotypic, dominant or recessive model in place of a multiplicative model, include the model option --genotypic, --dominant or --recessive, 
#respectively. 
#To specify interactions between covariates, and between SNPs and covariates, include the option --interaction.

# With genotypic the term is: Y = b0 + b1.ADD + b2.DOMDEV + b3.COV1 + b4.COV2 + e

# --genotypic will generate two extra tests per SNP, corresponding to: DOMDEV and GENO_2DF which represent a separate test of the dominance component 
# or a 2 df joint test of both additive and dominance. DOMDEV refers to a variable coded 0,1,0 for the three genotypes AA,Aa,aa, 
# i.e. representing the dominance deviation from additivity, rather specifying that a particular allele is dominant or recessive. 
#That is, the DOMDEV term is fitted jointly with the ADD term in a single model.
# Interpretation of results:
# That is, this represent coefficients from four terms in a multiple regression of disease on ADD, DOMDEV, COV1 and COV2 jointly. 
# The final test is a 2df test that tests the coefficients for ADD and DOMDEV together. 
# Importantly, the p-values for each line reflect the effect of the entity under the TEST column, not of the SNP whilst controlling 
# for that particular covariate. 

#vitd6
/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals \
	--pheno final_phenotype_data.tab --pheno-name vitd6 \
	--linear sex hide-covar --adjust --ci 0.95 --missing-phenotype -99999 \
	--covar final_phenotype_data.tab \
	--covar-name vitd0, incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, bmi0, calendar_age_ra, season_randomisation_2 \
	--out QCd_genotypes_vitd6_all_covar_vitd0 \
        --perm

#vitd12
/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals \
	--pheno final_phenotype_data.tab --pheno-name vitd12 \
	--linear sex hide-covar --adjust --ci 0.95 --missing-phenotype -99999 \
	--covar final_phenotype_data.tab \
	--covar-name vitd0, incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, bmi0, calendar_age_ra, season_randomisation_2 \
	--out QCd_genotypes_vitd12_all_covar_vitd0 \
        --perm

# vitd12 without correcting for baseline vitd (vitd0):
/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals \
        --pheno final_phenotype_data.tab --pheno-name vitd12 \
        --linear sex hide-covar --adjust --ci 0.95 --missing-phenotype -99999 \
        --covar final_phenotype_data.tab \
	--covar-name incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, bmi0, calendar_age_ra, season_randomisation_2 \
        --out QCd_genotypes_vitd12_all_covar_novitd0 \
        --perm


## Test with multiple covariates and interactions:

# TO DO: Test with --interaction? Unsupported in plink1.9? Flag doesnt' appear, use --gxe instead? (but limited functionality)
# Given both a quantitative phenotype and a case/control covariate (loaded with --covar) defining two groups, --gxe compares the regression 
# coefficient derived from considering only members of one group to the regression coefficient derived from considering only members of the other, 
# writing a report to plink.qassoc.gxe

# With --gxe to compare placebo vs treated. Note that covar index starts after FID and IID (so position 3 = 1):
# This command errors, see log file
/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals \
        --pheno final_phenotype_data.tab --pheno-name vitd12 \
        --assoc qt-means --gxe 56 --adjust --missing-phenotype -99999 \
        --covar final_phenotype_data.tab \
	--covar-name vitd0, incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, bmi0, calendar_age_ra, season_randomisation_2 \
        --out QCd_genotypes_vitd12_gxe_placebovsTx \
        --perm

# Linear regression test for vitd12 using only subgroups:
# All individuals:
/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals \
        --pheno final_phenotype_data.tab --pheno-name vitd12 \
        --linear sex hide-covar --ci 0.95 --adjust --missing-phenotype -99999 \
        --covar final_phenotype_data.tab \
	--covar-name vitd0, incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, bmi0, calendar_age_ra, season_randomisation_2 \
        --out QCd_genotypes_vitd12_all_samples_all_covar_vitd0 \
        --perm
# Get VD SNPs:
# TO DO: keep header for easier processing afterwards. Could also add column with source of SNP/test (ie filename?).
grep -wf VD_SNPs_GWAS_list.txt QCd_genotypes_vitd12_all_samples_all_covar_vitd0.assoc.linear.perm | sort -gk 12 \
        > vitd12_all_samples_all_covar_vitd0.assoc.linear.perm.results_SNPs_of_interest.txt


# Treated only:
/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals \
	--pheno final_phenotype_data.tab --pheno-name vitd12_2000_4000_only \
	--linear sex hide-covar --ci 0.95 --adjust --missing-phenotype -99999 \
	--covar final_phenotype_data.tab \
	--covar-name vitd0, incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, bmi0, calendar_age_ra, season_randomisation_2 \
	--out QCd_genotypes_vitd12_2000_4000_only_all_covar_vitd0 \
        --perm
# Get VD SNPs:
grep -wf VD_SNPs_GWAS_list.txt QCd_genotypes_vitd12_2000_4000_only_all_covar_vitd0.assoc.linear.perm | sort -gk 12 \
	> vitd12_2000_4000_only_all_covar_vitd0.assoc.linear.perm.results_SNPs_of_interest.txt

# 4000 only:
/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals \
        --pheno final_phenotype_data.tab --pheno-name vitd12_4000_only \
        --linear sex hide-covar --ci 0.95 --adjust --missing-phenotype -99999 \
        --covar final_phenotype_data.tab \
	--covar-name vitd0, incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, bmi0, calendar_age_ra, season_randomisation_2 \
        --out QCd_genotypes_vitd12_4000_only_all_covar_vitd0 \
        --perm
# Get VD SNPs:
grep -wf VD_SNPs_GWAS_list.txt QCd_genotypes_vitd12_4000_only_all_covar_vitd0.assoc.linear.perm | sort -gk 12 \
        > vitd12_4000_only_all_covar_vitd0.assoc.linear.perm.results_SNPs_of_interest.txt


# 2000 only:
/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals \
        --pheno final_phenotype_data.tab --pheno-name vitd12_2000_only \
        --linear sex hide-covar --ci 0.95 --adjust --missing-phenotype -99999 \
        --covar final_phenotype_data.tab \
	--covar-name vitd0, incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, bmi0, calendar_age_ra, season_randomisation_2 \
        --out QCd_genotypes_vitd12_2000_only_all_covar_vitd0 \
        --perm
# Get VD SNPs:
grep -wf VD_SNPs_GWAS_list.txt QCd_genotypes_vitd12_2000_only_all_covar_vitd0.assoc.linear.perm | sort -gk 12 \
        > vitd12_2000_only_all_covar_vitd0.assoc.linear.perm.results_SNPs_of_interest.txt

## Run the same for vitd6:
# All individuals:
/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals \
        --pheno final_phenotype_data.tab --pheno-name vitd6 \
        --linear sex hide-covar --ci 0.95 --adjust --missing-phenotype -99999 \
        --covar final_phenotype_data.tab \
	--covar-name vitd0, incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, bmi0, calendar_age_ra, season_randomisation_2 \
        --out QCd_genotypes_vitd6_all_samples_all_covar_vitd0 \
        --perm
# Get VD SNPs:
grep -wf VD_SNPs_GWAS_list.txt QCd_genotypes_vitd6_all_samples_all_covar_vitd0.assoc.linear.perm | sort -gk 12 \
        > vitd6_all_samples_all_covar_vitd0.assoc.linear.perm.results_SNPs_of_interest.txt


# Treated only:
/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals \
        --pheno final_phenotype_data.tab --pheno-name vitd6_2000_4000_only \
        --linear sex hide-covar --ci 0.95 --adjust --missing-phenotype -99999 \
        --covar final_phenotype_data.tab \
	--covar-name vitd0, incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, bmi0, calendar_age_ra, season_randomisation_2 \
        --out QCd_genotypes_vitd6_2000_4000_only_all_covar_vitd0 \
        --perm
# Get VD SNPs:
grep -wf VD_SNPs_GWAS_list.txt QCd_genotypes_vitd6_2000_4000_only_all_covar_vitd0.assoc.linear.perm | sort -gk 12 \
        > vitd6_2000_4000_only_all_covar_vitd0.assoc.linear.perm.results_SNPs_of_interest.txt

# 4000 only:
/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals \
        --pheno final_phenotype_data.tab --pheno-name vitd6_4000_only \
        --linear sex hide-covar --ci 0.95 --adjust --missing-phenotype -99999 \
        --covar final_phenotype_data.tab \
	--covar-name vitd0, incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, bmi0, calendar_age_ra, season_randomisation_2 \
        --out QCd_genotypes_vitd6_4000_only_all_covar_vitd0 \
        --perm
# Get VD SNPs:
grep -wf VD_SNPs_GWAS_list.txt QCd_genotypes_vitd6_4000_only_all_covar_vitd0.assoc.linear.perm | sort -gk 12 \
        > vitd6_4000_only_all_covar_vitd0.assoc.linear.perm.results_SNPs_of_interest.txt


# 2000 only:
/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals \
        --pheno final_phenotype_data.tab --pheno-name vitd6_2000_only \
        --linear sex hide-covar --ci 0.95 --adjust --missing-phenotype -99999 \
        --covar final_phenotype_data.tab \
	--covar-name vitd0, incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, bmi0, calendar_age_ra, season_randomisation_2 \
        --out QCd_genotypes_vitd6_2000_only_all_covar_vitd0 \
        --perm
# Get VD SNPs:
grep -wf VD_SNPs_GWAS_list.txt QCd_genotypes_vitd6_2000_only_all_covar_vitd0.assoc.linear.perm | sort -gk 12 \
        > vitd6_2000_only_all_covar_vitd0.assoc.linear.perm.results_SNPs_of_interest.txt


#vitd12 with interaction term all individuals:
/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals \
        --pheno final_phenotype_data.tab --pheno-name vitd12 \
        --linear interaction sex hide-covar --ci 0.95 --adjust --missing-phenotype -99999 \
	--covar final_phenotype_data.tab \
	--covar-name vitd0, incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, bmi0, calendar_age_ra, season_randomisation_2 \
        --out QCd_genotypes_vitd12_all_interaction_all_covar_vitd0 \
        --perm
grep -wf VD_SNPs_GWAS_list.txt QCd_genotypes_vitd12_all_interaction_all_covar_vitd0.assoc.linear.perm | sort -gk 12 \
        > vitd12_all_interaction_all_covar_vitd0.assoc.linear.perm.results_SNPs_of_interest.txt

# vitd12 with interaction term treated vs placebo:
/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals \
        --pheno final_phenotype_data.tab --pheno-name vitd12_2000_4000_only \
        --linear interaction sex hide-covar --ci 0.95 --adjust --missing-phenotype -99999 \
        --covar final_phenotype_data.tab \
	--covar-name vitd0, incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, bmi0, calendar_age_ra, season_randomisation_2 \
        --out QCd_genotypes_vitd12_2000_4000_only_interaction_all_covar_vitd0 \
        --perm
grep -wf VD_SNPs_GWAS_list.txt QCd_genotypes_vitd12_2000_4000_only_interaction_all_covar_vitd0.assoc.linear.perm | sort -gk 12 \
        > vitd12_2000_4000_only_interaction_all_covar_vitd0.assoc.linear.perm.results_SNPs_of_interest.txt

############

#Test without SNPs (i.e. linear regression without the genotypes):
# All individuals:
/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals \
        --pheno final_phenotype_data.tab --pheno-name vitd12 \
        --linear sex no-snp --ci 0.95 --missing-phenotype -99999 \
        --covar final_phenotype_data.tab \
	--covar-name vitd0, incident_fracture, incident_resp_infection, diabetes, heart_disease, copd_asthma, basemed_vitamind, currsmoker, bmi0, calendar_age_ra, season_randomisation_2 \
        --out QCd_genotypes_vitd12_all_samples_all_covar_vitd0_no_snps

# TO DO: command errors: With all-pheno test without SNPs for all phenotype data present (i.e. linear regression without the genotypes):
#/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals \
#        --pheno BEST-D_alternate_phenotypes_kit_ID_twice_recoded.txt --pheno-name vitd12 \
#        --linear sex no-snp --all-pheno --missing-phenotype -99999 \
#        --out QCd_genotypes_vitd12_no-snp_all_pheno

# all-pheno testing SNPs:
#/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals \
#	--pheno BEST-D_alternate_phenotypes_kit_ID_twice_recoded_Txsubgroups.txt \
#	--linear --all-pheno --adjust --ci 0.95 --missing-phenotype -99999 \
#	--out QCd_genotypes_vitd12_all_pheno

#Extract SNPs of interest from association analysis:

#/ifs/apps/apps/R-3.2.3/bin/Rscript extract_SNPs_assoc_test.R P140343-Results_FinalReport_clean_SNPs_autosome_individuals_vitd12_no_interaction.assoc.linear
# This R script doesn't keep duplicates (check match() function and re-run, easier with grep though:
#grep -wf VD_SNPs_GWAS_list.txt QCd_genotypes_vitd12_no_interaction.assoc.linear | sort -gk 12 \
#	> assoc_results_SNPs_of_interest.txt

#This produces two files (all matches, and matches with p-value cutoff) for each file tested.


## Extract subset of data for plotting/analysis in R:

# Convert to the plink ped (text) format if the file was binary:
/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals --recode --extract VD_SNPs_GWAS_list.txt \
        --out QCd_genotypes_VD_GWAS_SNPs

# Convert ped to csv with IDs in the first column and then genotypes:
# TO DO: modify so as to keep rsIDs as headers
cut -d " " -f2-2,7- --output-delimiter=, QCd_genotypes_VD_GWAS_SNPs.ped > QCd_genotypes_VD_GWAS_SNPs.csv

# Or convert to vcf format:
/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals --recode vcf --extract VD_SNPs_GWAS_list.txt \
	--out QCd_genotypes_VD_GWAS_SNPs

#########################

# Get errors, warnings, notes:
# Append to stdout: Move to ruffus and CGAT!
cat *.log | grep -i error
cat *.log | grep -i note
cat *.log | grep -i warning

#########################
#Data visualisation
#Quantile-quantile plot:
# Run this for all results. Modify scripts (move to Ruffus!):
/ifs/apps/apps/R-3.2.3/bin/Rscript plot-qqplot.R QCd_genotypes_vitd12_all_samples_all_covar_vitd0.assoc.linear.perm
/ifs/apps/apps/R-3.2.3/bin/Rscript plot-qqplot.R QCd_genotypes_vitd0_all_covar_linear.assoc.linear.perm

#Manhattan plot:
#Haploview -nogui -plink QCd_genotypes_vitd12_all_samples_all_covar_vitd0.assoc.linear.perm \
#	-map P140343-Results_FinalReport_clean_SNPs_autosome_individuals.bim \
#	-out QCd_genotypes_vitd12_all_samples_all_covar_vitd0

# Run genome-wide plots and regions/SNPs of interest (with all SNPs tested for the region):
/ifs/apps/apps/R-3.2.3/bin/Rscript plot-manhattan.R QCd_genotypes_vitd12_all_samples_all_covar_vitd0.assoc.linear.perm
/ifs/apps/apps/R-3.2.3/bin/Rscript plot-manhattan_SNPs_of_interest.R QCd_genotypes_vitd12_all_samples_all_covar_vitd0.assoc.linear.perm VD_SNPs_GWAS_list.txt
/ifs/apps/apps/R-3.2.3/bin/Rscript plot-manhattan_region_of_interest.R QCd_genotypes_vitd12_all_samples_all_covar_vitd0.assoc.linear.perm VD_SNPs_GWAS_list.txt 4


# Run as a candidate gene/region analysis:
# Look at LD relationships near potential hits 
# Get region +/- 250kb from peak SNP from PLINK ped file using --snp and --window:
# --recode HV gives Haploview format.
/ifs/apps/bio/plink-1.9329/plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals --snp rs2282679 \
	--window 250 --recode HV --out QCd_rs2282679_250kb

# Then load into Haploview using the linkage format

# Look at zoomed-in P-values for the region with LD values (for use in SNAP plot):
#plink2 --noweb --bfile P140343-Results_FinalReport_clean_SNPs_autosome_individuals \
#	--chr 4 --from-kb 71492666 --to-kb 71992666 \
#	--pheno BEST-D_alternate_phenotypes_kit_ID_twice_recoded_Txsubgroups.txt --pheno-name vitd12_2000_4000_only \
#	--linear sex hide-covar --ci 0.95 --adjust --missing-phenotype -99999  \
#        --covar BEST-D_alternate_phenotypes_kit_ID_twice_recoded_Txsubgroups.txt --covar-name bmi0,vitd0 \
#        --out QCd_genotypes_vitd12_2000_4000_only_all_covar_vitd0_chr4_71742666_250kb

#Change the header P to PValue 
#Go to SNAP plot website, choose "Plots" from upper right menu and plot a "Regional Association Plot"


 
## LD plot:


#########################
## Plotting with LocusZoom for regional associaiotn and LD plots:

# LocusZoom needs --metal format (two columns with specific headers, see:
# http://genome.sph.umich.edu/wiki/LocusZoom_Documentation
# http://genome.sph.umich.edu/wiki/LocusZoom_Standalone

# First convert Plink output (spaces/tabs get messed up, not sure why, use prettify tool or sed).
# Convert spaces to tabs, clean up plink output to make it easier to pass on:
#To perform the simplest reverse conversion (multiple spaces to one tab), you can use
#  cat [input filename] | tr -s ' ' '\t' > [output filename]
cat QCd_genotypes_vitd12_4000_only_all_covar_vitd0.assoc.linear.perm | tr -s ' ' '\t' \
	> QCd_genotypes_vitd12_4000_only_all_covar_vitd0.assoc.linear.perm.tsv

cat QCd_genotypes_vitd12_2000_only_all_covar_vitd0.assoc.linear.perm | tr -s ' ' '\t' \
        > QCd_genotypes_vitd12_2000_only_all_covar_vitd0.assoc.linear.perm.tsv

cat QCd_genotypes_vitd6_4000_only_all_covar_vitd0.assoc.linear.perm | tr -s ' ' '\t' \
        > QCd_genotypes_vitd6_4000_only_all_covar_vitd0.assoc.linear.perm.tsv

cat QCd_genotypes_vitd12_2000_only_all_covar_vitd0.assoc.linear.perm | tr -s ' ' '\t' \
        > QCd_genotypes_vitd6_2000_only_all_covar_vitd0.assoc.linear.perm.tsv


# This leaves a tab in the first column (as an empty column) though:
# To strip leading and trailing tabs and spaces, try:
#  cat [in] | sed 's/^[[:space:]]*//g' | sed 's/[[:space:]]*$//g' > [out]
#cat QCd_genotypes_vitd12_2000_4000_only_all_covar_vitd0.assoc.linear.perm.tsv | sed 's/^[[:space:]]*//g' | \
#	sed 's/[[:space:]]*$//g' > QCd_genotypes_vitd12_2000_4000_only_all_covar_vitd0.assoc.linear.perm.tsv
# TO DO: didn't work, deletes content. Left with empty first column for now as doesn't affect.

#########
# TO DO: clean up and automate, pass top hits from each association result to locusZoom:

# Plot with locusZoom:
#cat QCd_genotypes_vitd12_2000_4000_only_all_covar_vitd0.assoc.linear.perm.tsv | cut -f3,13 | locuszoom --metal - \
#	--markercol SNP --pvalcol P --refsnp rs2282679 --flank 500kb \
#	--pop EUR --build hg19 --source 1000G_March2012 \
#	--gwas-cat whole-cat_significant-only --build hg19 \
#	--prefix QCd_genotypes_vitd12_2000_4000_only_all_covar_vitd0.assoc.linear.perm

# Run with a batch file to plot many SNPs/regions with --hitspec VD_SNPs_GWAS_LD_locusZoom_batch.txt
# Errors with I/O operation, send bug report...

#cat QCd_genotypes_vitd12_2000_4000_only_all_covar_vitd0.assoc.linear.perm.tsv | cut -f3,13 | locuszoom --metal - \
#        --markercol SNP --pvalcol P --hitspec VD_SNPs_GWAS_LD_locusZoom_batch.txt \
#        --pop EUR --build hg19 --source 1000G_March2012 \
#        --gwas-cat whole-cat_significant-only --build hg19 \
#        --prefix QCd_genotypes_vitd12_2000_4000_only_all_covar_vitd0.assoc.linear.perm


# Run manually for top SNPs:
# vitd6 and vitd12 at 2000 only adjusting for covariates:
cat QCd_genotypes_vitd12_2000_only_all_covar_vitd0.assoc.linear.perm.tsv | cut -f3,13 | /ifs/apps/bio/locusZoom/locuszoom --metal - \
        --markercol SNP --pvalcol P --refsnp rs222047 --flank 500kb \
        --pop EUR --build hg19 --source 1000G_March2012 \
        --gwas-cat whole-cat_significant-only --build hg19 \
        --prefix QCd_genotypes_vitd12_2000_only_all_covar_vitd0.assoc.linear.perm

cat QCd_genotypes_vitd6_2000_only_all_covar_vitd0.assoc.linear.perm.tsv | cut -f3,13 | /ifs/apps/bio/locusZoom/locuszoom --metal - \
        --markercol SNP --pvalcol P --refsnp rs222047 --flank 500kb \
        --pop EUR --build hg19 --source 1000G_March2012 \
        --gwas-cat whole-cat_significant-only --build hg19 \
        --prefix QCd_genotypes_vitd6_2000_only_all_covar_vitd0.assoc.linear.perm

# vitd12 2000 only p-value ~0.08, vitd6 2000 only significant:
cat QCd_genotypes_vitd6_2000_only_all_covar_vitd0.assoc.linear.perm.tsv | cut -f3,13 | /ifs/apps/bio/locusZoom/locuszoom --metal - \
        --markercol SNP --pvalcol P --refsnp rs7041 --flank 500kb \
        --pop EUR --build hg19 --source 1000G_March2012 \
        --gwas-cat whole-cat_significant-only --build hg19 \
        --prefix QCd_genotypes_vitd12_2000_4000_only_all_covar_vitd0.assoc.linear.perm

# vitd12 2000 only p-value ~0.08, vitd6 2000 only significant:
cat QCd_genotypes_vitd6_2000_only_all_covar_vitd0.assoc.linear.perm.tsv | cut -f3,13 | /ifs/apps/bio/locusZoom/locuszoom --metal - \
        --markercol SNP --pvalcol P --refsnp rs4588 --flank 500kb \
        --pop EUR --build hg19 --source 1000G_March2012 \
        --gwas-cat whole-cat_significant-only --build hg19 \
        --prefix QCd_genotypes_vitd6_2000_only_all_covar_vitd0.assoc.linear.perm

# vitd6 4000 only, significant:
cat QCd_genotypes_vitd6_4000_only_all_covar_vitd0.assoc.linear.perm.tsv | cut -f3,13 | /ifs/apps/bio/locusZoom/locuszoom --metal - \
        --markercol SNP --pvalcol P --refsnp rs2282621 --flank 500kb \
        --pop EUR --build hg19 --source 1000G_March2012 \
        --gwas-cat whole-cat_significant-only --build hg19 \
        --prefix QCd_genotypes_vitd6_4000_only_all_covar_vitd0.assoc.linear.perm

cat QCd_genotypes_vitd6_4000_only_all_covar_vitd0.assoc.linear.perm.tsv | cut -f3,13 | /ifs/apps/bio/locusZoom/locuszoom --metal - \
        --markercol SNP --pvalcol P --refsnp rs11023350 --flank 500kb \
        --pop EUR --build hg19 --source 1000G_March2012 \
        --gwas-cat whole-cat_significant-only --build hg19 \
        --prefix QCd_genotypes_vitd6_4000_only_all_covar_vitd0.assoc.linear.perm


# TO DO: Add a file with markers to highlight: --denote-markers-file VD_SNPs_GWAS_list.txt 
# This file should be tab delimited and have the headers:
#snp 	string 	color
#rs231362 	GWAS 	blue
#rs163184 	Conditional 	purple

# TO DO: If a list of fine-mapped SNPs is available for the region (a 'credible set'), 
# these can be added as a plotting
# option as ' fineMap="my_finemapping_results.txt" '. See http://genome.sph.umich.edu/wiki/LocusZoom_Standalone
# Check: http://www.nature.com/ng/journal/v44/n12/full/ng.2435.html

# TO DO: add a bed track such as ENCODE regions: --bed-tracks <your bed file> 
# The BED file should have at least 4 columns: the first 3 for chr/start/end, and the 4th column for the label
# Colours can be added, follow the full BED format.

# Other options, customisation, etc: http://genome.sph.umich.edu/wiki/LocusZoom#Commonly_Used_LocusZoom_Options 


#########################
#Multiple testing adjustment

#########################

#########################
#Correction for population stratification


#########################


#########################

#TO DO: PCA as in eQTLs for VD levels, CRP, etc.
#TO DO: Assoc. analysis with VitD supp. change in levels as phenotype.

#########################


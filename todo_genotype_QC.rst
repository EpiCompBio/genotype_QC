################################
Plan for genotype QC for Airwave
################################

:Date: 22 Dec 2016
:Authors: Gao He, Vangelis, Daivd M, Antonio B


Current state:
##############

- The plan is to create a QC pipeline based on Plink and other tools that can be automated using Ruffus and CGAT tools. 

- Based on Gao's and other work the genotype data looks very good in general. 

- CGAT, ruffus, drmaa are in place and generally working with a few bumps.

- We have most of the scripts we need although they require a lot of cleaning up.

- Gao, Vangelis and I have discussed and agreed the steps to follow. These basically mirror de UK Biobank set up.

- There is a pipeline in progress in the src folder with a dummy report in the report folder (based on CGAT's quickstart pipeline).


Future work:
############

- After the QC steps we will discuss the need for an imputation pipeline and what to do with the different platforms. i.e. Merge them, run as validation, as meta-analysis, etc. This may be project/question dependent though.

- For imputation it would be ideal to do it ourselves. This needs person programming time and accessing the latest reference panel (we should try to get the HRC, Vangelis is looking into this).



TO DO main tasks:
#################

- Most of the scripts (from Gao and I and third party) are ready but need cleaning up, making sure they are callable and put into a pipeline.

- I've marked what needs to be done below with "TO DO" and a short explanantion, mostly either creating the stand-alone script for that step or checking the existing one and transferring the command into the ruffus pipeline or cleaning it up/turning into a script.


Notes for Airwave data: 
#######################

.. todo::
::

TO DO: Location of basic phenotype data
TO DO: location of current genotype data
TO DO: location of WTCHG Core Facility reports

.. note:: There are three platforms: Exome, CoreExome and an Affy chip

.. todo:: 
	TO DO: location of each

.. note:: Plates 1-27 are batch 1, 28-53 are batch 2.


Notes from Gao's work:
######################

Gao has done cleaning up, formatting, sanity checks, etc. already. 
See emails

-----

Subject: RE: airwave QC discussion
From: "Gao, He" <h.gao@imperial.ac.uk>
Date: 19/09/2016, 17:17
To: "Berlanga, Antonio J" <a.berlanga@imperial.ac.uk>
CC: "Evangelou, Evangelos" <e.evangelou@imperial.ac.uk>

- updated the alleles (from A/B to TOP/BOT annotation), strand and genome build. 
- updated new plink files
- calculated the sample and SNP missingness and HWE p-values PER BATCH. Sample missingness is very low in general. SNP missingness looks ok for all batches except for n23 and n46. (I think when calculating the SNP call rates, we need to exclude the non-autosomal SNPs.)
 
TO DO: move these scripts here:
- The script (2_alleleBuildUpdate.pbs, 3_quickCheck.pbs) and files are all in the folder /groupvol/med-bio/epiUKB/Airwave/coreExome_genotype/plinkFiles. They are generated for information purpose and can be removed if we don’t need them later.
 
- fixed the four batches and the plink files have been generated. In case you need to run any loop there should be no problem now.
 
- calculated the number of samples and some full batches do not have 288 samples in the intensity files (and therefore also in the plink files).

.. todo::
	TO DO: more emails missing, update/transfer info here for clarity and follow-up

-----

test data: use first few batches for example, see bed, bim and fam files in:
/groupvol/med-bio/epiUKB/Airwave/coreExome_genotype/plinkFiles/


PIPELINE PLAN
#############

Files to check to pull out commands for most of the steps below:
|	plink_preprocessing.txt
|	QC_plink_individuals.sh
|	markers_QC_Airwave.sh

These scripts were run as QC of markers first, then individual samples. Steps in this pipeline follow the UK Biobank protocol (which goes back and forth between markers and individuals as it first creates a homogeneous group [based on ethnicity] from where to pull out high quality genetic markers which are not confounded by population stratification.

.. note::
	Keep scripts, modules, pipelines separate.

.. todo::
	TO DO: scan/ppt pipeline workflow plus notes

-----

A. Pre-QC steps, GenomeStudio to plink, hg19 liftover, flip strand:

	TO DO: load into pipeline by calling each script or function. Needs a if/else decision (if illumina, convert to xxx, if affy do xxx, else error):

	1. GenomeStudio to plink: by zcall script:
		Script: /groupvol/med-bio/epiUKB/Airwave/coreExome_zcall/zcall_v3.4/convertReportToTPED.py
		Job submission script: /groupvol/med-bio/epiUKB/Airwave/coreExome_genotype/plinkFiles/1_convertReportToTPED.pbs
		Result files: /groupvol/med-bio/epiUKB/Airwave/coreExome_genotype/plinkFiles

	2. Convert from AB allele to illumina TOP/BOT annotation: by plink, using Wrayner's annotation files
		Strand files: /groupvol/med-bio/epiUKB/Airwave/strandFiles
		(from http://www.well.ox.ac.uk/~wrayner/strand/)
		Command: plink --noweb --bfile --update-alleles humancoreexome-12v1-1_a.update_alleles.txt --make-bed --out

	3. Update genome build: hg19/build 37 liftover: by plink, using Wrayner's annotation files, also handles strand
		This includes updating a few attributes (chromosome, position, strand flipping etc)
		Script: http://www.well.ox.ac.uk/~wrayner/strand/update_build.sh

-----

B. Allele frequency report with proportions:
	TO DO write commands into ruffus pipeline, e.g. (see also sh scripts above):
	plink2 --bifle xxx --freq
	cat plink.frq | tr -s ' ' '\t' | cut -f 4 | grep A | wc -l # First column is a tab, so fourth is A1

-----

#. Select homogeneous set of samples to use as set for marker QC (PCA based, with automatic selection using e.g. 'aberrant' R package. This is to avoid artefacts from population structure. Excluded samples are later re-introduced.):
	http://bioinformatics.oxfordjournals.org/content/28/1/134.full.pdf+html
	Use summary statistics, and/or: missingness, ancestry, probe intensity, gender separately:
	TO DO write commands into ruffus pipeline:
		- Merge plates first
		
	TO DO write commands into ruffus pipeline (see scripts above although PCA tool needs changing to FlashPCA probably as older tools won't run on large number of samples):
		- Run PCA against 1000G (or Hapmap) as in UKB appendix 1 (requires using plink MAF >5%, HWE 10^-6, etc for Hapmap or 1000G, then projecting onto these)
		
	TO DO write script to wrap aberrant and make it callable from CLI within pipeline:	
		- aberrant with lambda set to 20 for ancestry PC1 and PC2 as summary stats

-----

#. Per batch marker QC (plink commands; drop failing SNPs from all plates):
	- TO DO write script for this, needs loop calling batch 1 vs all other batches, then batch 2 vs all other batches, etc. with parameters (eg p-values and all the criteria below) can be set by user:
		+ Exclude monomorphic SNPs
		+ Genotype call rate (<98%)
		+ Genotype frequency consistency across batches (Fisher's exact test p-value <10^-12)
		+ Allele frequency consistency versus reference panel (eg Hapmap, Fisher's exact test p-value <10^-12)
		+ Hardy Weinberg equilibrium (p-value <10^-12)

-----

#. Plate/batch PCA (visual outlier detection check)
	TO DO clean up commands from above and plotting script for this (may need substantial re-working with tools that take thousands of samples, check notes/discuss)

-----

#. Plate/batch merge
	TO DO write scripts/commands

-----

#. Visual test of genotype calls in cluster plots (bin by MAF, pick random subset)
	TO DO write scripts for this: Gao has plotted these before and I think has scripts. Obviously can't check thousands of SNPs visually svo either use a random pick (e.g. grab 20 or whatever is plottable) or better grab top 10 highest quality SNPs, bottom 10, 10 failed SNPs, 10 at MAF > 10%, 10 at 1-5%, 10 <1%, etc. The aim is to have some visual sanity check of the raw data for some of the markers.

-----

#. Pooled sample QC (all samples; based on high quality set of markers from above; plink commands):
	TO DO these are plink commands that can be put directly into the ruffus pipeline with a PARAMS config option so user can set different cut-offs (these PARAMS and config file are standard for CGAT pipelines):
     - Run with autosomal SNPs only
     - Heterozygosity (standard deviation > +/- 3) and genotype failure rates per individual (>5%)
     - Relatedness between individuals (IBD cut-off >0.185)
     - Gender mis-identification check

-----


#. VCF check sanity (strand, problematic SNPs, etc.)
TO DO look up tools and insert command into Ruffus, these already exist, plink2 has commands for this.



Future work
###########

The output of the genotpe QC pipeline should be input for:

- Imputation + post-imputation QC (discuss this at a later stage)

- If several platforms, re-run whole script with each platform as separate batch

- Final set of QC'd and imputed individuals and markers can then be processed depending on specific question, eg:
     + Filter out MAF <1%
     + Ancestry PCA (e.g. vs Hapmap populations), filter out on PC2 values <0.072
     + Filter out mitochondrial and sex chromosomes

- VCF check sanity (strand, problematic SNPs, etc.)


References
##########

General protocols and references:
http://www.ukbiobank.ac.uk/wp-content/uploads/2014/04/UKBiobank_genotyping_QC_documentation-web.pdf
http://www.nature.com/nprot/journal/v5/n9/pdf/nprot.2010.116.pdf
http://www.nature.com/nprot/journal/v10/n9/pdf/nprot.2015.077.pdf
http://www.nature.com/ng/journal/vaop/ncurrent/pdf/ng.3656.pdf


Also see:
Quality control and conduct of genome-wide association meta-analyses
http://www.nature.com/nprot/journal/v9/n5/full/nprot.2014.071.html

Basic statistical analysis in genetic case-control studies
http://www.nature.com/nprot/journal/v6/n2/abs/nprot.2010.182.html

Further references in:
https://github.com/EpiCompBio/genotype_tools/blob/master/src/pipeline_genotype_QC.py


Downstream annotation
#####################

.. todo:: 
	move this to the next pipeline

DEPICT Biological interpretation of genome-wide association studies using predicted gene functions.
http://www.ncbi.nlm.nih.gov/pubmed/25597830?dopt=Abstract&holding=npg

# Hello

##############################################################################
#
#   Imperial College London
#
#   $Id$
#
#   Copyright (C) 2017 Antonio Berlanga-Taylor
#
#   This program is free software; you can redistribute it and/or
#   modify it under the terms of the GNU General Public License
#   as published by the Free Software Foundation; either version 3
#   of the License, or (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program; if not, write to the Free Software
#   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
###############################################################################

#########################
#Genotype QC pipeline
#Pre-processing of raw genotype data for use as plink input and QC to eliminate markers and individuals.
#Antonio J Berlanga-Taylor
#########################


"""===========================
Pipeline genotype QC
===========================

:Author: Antonio Berlanga-Taylor
:Release: $Id$
:Date: |today|
:Tags: Python

.. Replace the documentation below with your own description of the
   pipeline's purpose

Overview
========

##Quality control of subject and SNP genotyping data
# Note this is now using plink 1.9, which is currently in beta

This pipeline pre-processes genome-wide genotype data from a microarray (Illumina, Affymetrix) and carries out quality control
processing. It outputs varies plots, tables and a final QC'd file for downstream analysis (GWAS, imputation, etc.).

See the files :file:``pipeline.ini` and :file:`conf.py` to configure with different parameters and reporting.

The pipeline is based on the UK Biobank protocol and other methods. See references. 

It requires CGAT tools (pipelines and scripts).

See notes and further information in:
https://github.com/EpiCompBio/genotype_tools/blob/master/todo_genotype_QC.rst

Anderson et al. 2010 protocol: http://www.nature.com/nprot/journal/v5/n9/pdf/nprot.2010.116.pdf
Winkler et al. 2014 protocol (meta-analysis of GWAS): http://www.nature.com/nprot/journal/v9/n5/pdf/nprot.2014.071.pdf
Plink tutorial: http://pngu.mgh.harvard.edu/~purcell/plink/tutorial.shtml

Usage
=====

See :ref:`PipelineSettingUp` and :ref:`PipelineRunning` on general
information how to use CGAT pipelines.

Configuration
-------------

The pipeline requires a configured :file:`pipeline.ini` file.
CGATReport report requires a :file:`conf.py` and optionally a
:file:`cgatreport.ini` file (see :ref:`PipelineReporting`).

Default configuration files can be generated by executing:

   python <srcdir>/pipeline_genotype_QC.py config

Input files
-----------

Genotype file as provided by Illumina, Affymetrix or converted to Plink formats.

TO DO: For pre-processing input is an Illumina SNP genotype file

Requirements
------------

The pipeline requires various tools to be installed:

Requirements:

Default CGAT setup: 
* CGATPipelines
* cgat tools

* Plink 1.90
* R
* FlashPCA: ftp://pricelab:pricelab@ftp.broadinstitute.org/EIGENSOFT/EIG6.0.1.tar.gz
* aberrant: http://www.well.ox.ac.uk/~spencer/Aberrant/aberrant-manual.pdf
* Some plotting scripts and others from Anderson et al 2010 protocol

TODO: add versions

Pipeline output
===============

Outputs a genetic marker and individual QC'd file in plink's format plus various descriptive plots and tables 
in a simple report.


Glossary
========

.. glossary::


Code
====

"""

#####################################################################
#####################################################################
# Several of the following functions are from CGATPipelines
#####################################################################

from ruffus import *

import sys
import os
import sqlite3

import pandas as pd
import numpy as np

import CGAT.Experiment as E
import CGATPipelines.Pipeline as P
import CGAT.Database as DB

# load options from the config file:

PARAMS = P.getParameters(
    ["%s/pipeline.ini" % os.path.splitext(__file__)[0],
     "../pipeline.ini",
     "pipeline.ini"])

# TO DO: check the ini file and configure appropriately, check annot pipeline:

# add configuration values from associated pipelines
#
# 1. pipeline_annotations: any parameters will be added with the
#    prefix "annotations_". The interface will be updated with
#    "annotations_dir" to point to the absolute path names.

PARAMS.update(P.peekParameters(
    PARAMS["annotations_dir"],
    "pipeline_annotations.py",
    on_error_raise=__name__ == "__main__",
    prefix="annotations_",
    update_interface=True))


# if necessary, update the PARAMS dictionary in any modules file.
# e.g.:
#
# import CGATPipelines.PipelineGeneset as PipelineGeneset
# PipelineGeneset.PARAMS = PARAMS
#
# Note that this is a hack and deprecated, better pass all
# parameters that are needed by a function explicitely.


#####################################################################
#####################################################################
# Check file inputs, error if none or not the correct file formats
#####################################################################

# TO DO: check PARAMS calling and pipeline.ini file, change this for os.sys check, not params:

if PARAMS["input"].lower() == "bed":
    suffix_pattern = "*.bed"
else:
    raise ValueError('''Binary plink input files are needed (bed, bim and fam, see 
    https://www.cog-genomics.org/plink2/formats#bed ; If you only have a .bed generate a .bim 
    and a dummy .fam''')

# Check if there is a bim and fam files as well:
# TO DO: How to handle downstream if there is only a bed file (and no bim and fam)? Error here for now, 
# easier to leave to user.

#else:
#raise ValueError("No bim and/or fam files detected. Check suffixes are OK, If you only have a .bed file generate a .bim 
#    from this one and a dummy .fam in order to run this pipeline.")

# Check there is at least one input file (one bed, there should be bim and fam as well though):
if len(SAMPLE_FILES) == 0:
    raise ValueError("No input files in the working directory")

########################
########################
# CGATPipeline function:
########################

# -----------------------------------------------
# Utility functions

def connect():
	'''
    utility function to connect to database.

    Use this method to connect to the pipeline database.
    Additional databases can be attached here as well.

    Returns an sqlite3 database handle.
	'''

    dbh = sqlite3.connect(PARAMS["database_name"])
    statement = '''ATTACH DATABASE '%s' as annotations''' % (
        PARAMS["annotations_database"])
    cc = dbh.cursor()
    cc.execute(statement)
    cc.close()

    return dbh


#########################################
#########################################
# Set up some general job mission options
#########################################

# Set-up per job, leave global definitions based on CGAT's priority handling, in beta though but current form works:
#https://github.com/CGATOxford/CGATPipelines/pull/254

#to_cluster = True
#job_threads = 
#job_memory = 



#########################
#########################
# Specific pipeline tasks
#########################

'''
General pipeline steps:
-----

A. Pre-QC steps, GenomeStudio to plink, hg19 liftover, flip strand:

# TO DO: load into pipeline by calling each script or function. Needs a if/else decision (if illumina, convert to xxx, if affy do xxx, else error):

#Data management see:
#Convert files from raw genotypes to Plink's lgen, map and ped formats:
#http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#ped
#https://www.biostars.org/p/10332/
#http://pngu.mgh.harvard.edu/~purcell/plink/data.shtml#map

From Gao, check if we can run the simple commands I have instead (to reduce external script dependency):

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
'''

# TO DO: check call to params, active if illumina = true:
@active_if(PARAMS["illumina"])
@transform(("*.txt",
	    suffix(".txt"),
	    "")
def convertIllumina(infile, outfile):
    '''
    Convert Illumina's Beadstudio output (e.g. 'xxx_FinalReport.txt') into 
    Plink's lgen, map, ped and fam formats.
    Check Plink output for results: log file, nof (nofounder) report, nosex (ambiguous sex) report.
    '''

	# Convert from raw to lgen
	# Convert from raw to map
	# Use lgen and map to obtain fam (lgen file has no paternal/maternal IDs, sex or phenotype information 
	# so are set to missing, see below for generating new fam file and alternate phenotype file)
	# Use lgen, fam and map to generate ped

    # the command line statement we want to execute:

    statement = '''
        cat %(infile)s | cut -f1,2,17,18 | sed "1,10d" | sed "s/-/0/g" | awk '{print $2,$2,$1,$3,$4}' > %(outfile)s.lgen;
        checkpoint;
        cat %(infile)s | cut -f1,19,20 | sed "1,10d" | awk '{print $2,$1,"0",$3}' | sort -u > %(outfile)s.map;
        checkpoint;
        perl -ane '{print "$F[0] $F[0] 0 0 0 -9\n"}' %(outfile)s.lgen | sort -u -k1,1 > %(outfile)s.fam;
        checkpoint;
        plink2 --noweb --lfile %(outfile)s --recode --out %(outfile)s.ped;
        checkpoint;
        touch %(outfile)s;
        '''

	# Generate binary files with plink2 (faster reading):
	# To disable the automatic setting of the phenotype to missing if the individual has an ambiguous 
	# sex code, add the --allow-no-sex option.
	# Sanity check with plink2 to see if file is intact and generate some summary stats. 
	# Results should be the same as log file generated from converting to binary file
    statement = '''
        plink2 --noweb --file %(outfile)s --make-bed --out %(outfile)s;
        checkpoint;
        plink2 --noweb --bfile %(outfile)s;
        '''

    # execute command in variable statement.
    # The command will be sent to the cluster.  The statement will be
    # interpolated with any options that are defined in in the
    # configuration files or variable that are declared in the calling
    # function.  For example, %(infile)s will we substituted with the
    # contents of the variable "infile".
    P.run()

'''
# Dummy function, loads data to a database, this is a CGAT function:
@transform(countWords,
           suffix(".counts"),
           "_counts.load")
def loadWordCounts(infile, outfile):
    " '''load results of word counting into database.''' " 
    P.load(infile, outfile, "--add-index=word")

# TO DO: Also see (from Steve S. pipeline):
@merge(counts,
       "counts.dir/counts.load")
def loadCounts(infiles, outfile):

    P.concatenateAndLoad(infiles, outfile,
                         regex_filename=".*/(.*).counts.gz",
                         has_titles=False,
                         cat="track",
                         header="track,ID,counts",
                         options='-i "ID"',
                         job_memory=PARAMS["xxx"])


# TO DO: Also see (from Steve S. pipeline):

@transform(someResults,
           suffix(".log"),
           ".load")
def loadSuperResult(infile, outfile):
    '''load the xxx results/processed data table from xxx function into the database'''

    my_table = os.path.dirname(infile) + "/my_name.my_table"

    P.load(my_table, outfile, 
	   options='-i "tracking_id"')

'''
           
@follows(xxx)
def preProcessIllumina():
    '''preProcessIllumina target'''
    pass


#########################		   

'''
-----

B. Allele frequency report with proportions:
	TO DO write commands into ruffus pipeline, e.g. (see also sh scripts above):
	plink2 --bifle xxx --freq
	cat plink.frq | tr -s ' ' '\t' | cut -f 4 | grep A | wc -l # First column is a tab, so fourth is A1
'''

# Get frequency stats:
#plink --bfile chr22_Airwave_CPMG_Plasma --freq --out test_chr22.test

@follows(xxx)
def alleleFreq():
    '''alleleFreq target'''
    pass


'''
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

'''

# Run loop to convert vcf to binary and subset individuals:
# Bash loop to process with qsub (use #Ruffus!):
#qsub PBS_plink_subset_loop.sh      
           
# Create list of subsetted files:
#ls -1 *.subset.bed > list_subsetted_files.txt

# Merge all resulting plink binaries into one file:
#plink --make-bed --bfile [?] --merge-list list_subsetted_files.txt --out Airwave_imputed_metabolomics_subset



           
@follows(xxx)
def homogeneousSet():
    '''homogeneousSet target'''
    pass


'''
-----

#. Per batch marker QC (plink commands; drop failing SNPs from all plates):
	- TO DO write script for this, needs loop calling batch 1 vs all other batches, then batch 2 vs all other batches, etc. with parameters (eg p-values and all the criteria below) can be set by user:
		+ Exclude monomorphic SNPs
		+ Genotype call rate (<98%)
		+ Genotype frequency consistency across batches (Fisher's exact test p-value <10^-12)
		+ Allele frequency consistency versus reference panel (eg Hapmap, Fisher's exact test p-value <10^-12)
		+ Hardy Weinberg equilibrium (p-value <10^-12)

'''


@follows(xxx)
def markerQC():
    '''markerQC target'''
    pass

'''
-----

#. Plate/batch PCA (visual outlier detection check)
	TO DO clean up commands from above and plotting script for this (may need substantial re-working with tools that take thousands of samples, check notes/discuss)

'''

@follows(xxx)
def PCA():
    '''PCA target'''
    pass

'''
-----

#. Plate/batch merge
	TO DO write scripts/commands

'''

@follows(xxx)
def mergePlates():
    '''mergePlates target'''
    pass

'''
-----

#. Visual test of genotype calls in cluster plots (bin by MAF, pick random subset)
	TO DO write scripts for this: Gao has plotted these before and I think has scripts. Obviously can't check thousands of SNPs visually svo either use a random pick (e.g. grab 20 or whatever is plottable) or better grab top 10 highest quality SNPs, bottom 10, 10 failed SNPs, 10 at MAF > 10%, 10 at 1-5%, 10 <1%, etc. The aim is to have some visual sanity check of the raw data for some of the markers.

'''



@follows(xxx)
def SNPcluster():
    '''SNPcluster target'''
    pass

'''
-----

#. Pooled sample QC (all samples; based on high quality set of markers from above; plink commands):
	TO DO these are plink commands that can be put directly into the ruffus pipeline with a PARAMS config option so user can set different cut-offs (these PARAMS and config file are standard for CGAT pipelines):
     - Run with autosomal SNPs only
     - Heterozygosity (standard deviation > +/- 3) and genotype failure rates per individual (>5%)
     - Relatedness between individuals (IBD cut-off >0.185)
     - Gender mis-identification check

'''


@follows(xxx)
def sampleQC():
    '''sampleQCsanity target'''
    pass

'''
-----

#. VCF check sanity (strand, problematic SNPs, etc.)
TO DO look up tools and insert command into Ruffus, these already exist, plink2 has commands for this.

'''

@follows(xxx)
def VCFsanity():
    '''VCFsanity target'''
    pass

##################################################################
# ---------------------------------------------------
# Generic pipeline tasks
@follows(xxx)
def full():
    pass

##################################################################
# ---------------------------------------------------
# Generic pipeline tasks for CGATReport:

@follows(mkdir("report"))
def build_report():
    '''build report from scratch. Any existing report will be overwritten.'''

    E.info("starting report build process from scratch")
    P.run_report(clean=True)


@follows(mkdir("report"))
def update_report():
    '''update report.

    This will update a report with any changes inside the report
    document or code. Note that updates to the data will not cause
    relevant sections to be updated. Use the cgatreport-clean utility
    first.
    '''

    E.info("updating report")
    P.run_report(clean=False)


@follows(update_report)
def publish_report():
    '''publish report in the CGAT downloads directory.'''

    E.info("publishing report")
    P.publish_report()

##################################################################
# ---------------------------------------------------
# See Steve's way of reporting:
# https://github.com/snsansom/scseq/blob/master/pipelines/pipeline_scrnaseq.py

# --------------------- < generic pipeline tasks > -------------------------- #

'''
@follows(mkdir("notebook.dir"))
@transform(glob.glob(os.path.join(os.path.dirname(__file__),
                                  "pipeline_notebooks",
                                  os.path.basename(__file__)[:-len(".py")],
                                  "*")),
           regex(r".*/(.*)"),
           r"notebook.dir/\1")
def notebooks(infile, outfile):
    ' '''Utility function to copy the notebooks from the source directory
       to the working directory''' '

    shutil.copy(infile, outfile)


@follows(quantitation, qc, notebooks)
def full():
    pass

print sys.argv

if __name__ == "__main__":
	sys.exit(P.main(sys.argv))

# See:
# https://github.com/snsansom/scseq/tree/master/pipelines/pipeline_notebooks/pipeline_scrnaseq
'''

##################################################################
# ---------------------------------------------------
# The end:

print sys.argv

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))

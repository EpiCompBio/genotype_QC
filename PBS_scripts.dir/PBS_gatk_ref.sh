#!/bin/sh

## This follows GATK's instructions to format a reference for use in other GATK commands. See:
##https://www.broadinstitute.org/gatk/guide/article?id=2798

## Input is a reference fasta file:
## Download from (April 2016):
##http://hgdownload.cse.ucsc.edu/goldenPath/hg38/bigZips/

## Outputs are various bwa files, samtools index and picard's dictionary.

## PBS job submission script:

#Job parameters:
#PBS -l walltime=00:10:00
#PBS -l select=1:ncpus=16:mem=8gb
##PBS -l select=NN:ncpus=MM:tmpspace=500gb
## This tells the batch manager to use NN nodes with MM cpus and with at least 500gb of free space on the temporary directory.

#Use current working directory
##$ -cwd

#Save standard error and out to files:
#$ -e stderr.file
#$ -o stdout.file

#Module commands:
module load gatk/3.5
module load bio-bwa
module load samtools/1.2
module load picard/2.1.1

#INFILE=hg38.fa
INFILE=chr22.fa
#OUTFILE=hg38.dict
OUTFILE=chr22.dict

SLOTS=16
#CONFIG=.ncbirc

#File management:
cp $PBS_O_WORKDIR/$INFILE $TMPDIR
##cp $PBS_O_WORKDIR/$CONFIG $TMPDIR

#Command:
bwa index -a bwtsw $INFILE
gunzip $INFILE
samtools faidx $INFILE
picard CreateSequenceDictionary REFERENCE=$INFILE OUTPUT=$OUTFILE 

#Notify by email:
##PBS -m abe


#File management:
cp $INFILE.* $PBS_O_WORKDIR
cp $OUTFILE $PBS_O_WORKDIR

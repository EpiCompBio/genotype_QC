#!/bin/sh

#Job parameters:
#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=16:mem=8gb
#PBS -q med-bio

#Use current working directory
##$ -cwd

#Save standard error and out to files:
##$ -e stderr.file
##$ -o stdout.file

#Module commands:
module load bcftools/1.2

INFILE=chr22_header2_head10000.dose.vcf
##chr22_header.dose.vcf.gz
OUTFILE=chr22_Airwave_CPMG_Plasma.vcf
SAMPLE_IDs=Airwave_CPMG_Plasma.txt_sample_names_FID_IID.txt
SLOTS=16
##CONFIG=.ncbirc

#File management:
cp $PBS_O_WORKDIR/$INFILE $TMPDIR
cp $PBS_O_WORKDIR/$SAMPLE_IDs $TMPDIR
##cp $PBS_O_WORKDIR/$CONFIG $TMPDIR

#Command:
bcftools view --samples-file $SAMPLE_IDs --force-samples $INFILE -O v -o $OUTFILE

#Notify by email:
##PBS -m abe


#File management:
cp $OUTFILE $PBS_O_WORKDIR

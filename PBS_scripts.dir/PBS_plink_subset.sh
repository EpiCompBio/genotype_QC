#!/bin/sh

#Job parameters:
#PBS -l walltime=05:00:00
#PBS -l select=1:ncpus=16:mem=10gb
#PBS -q med-bio 

#Use current working directory
##$ -cwd

#Save standard error and out to files:
##$ -e stderr.file
##$ -o stdout.file

#Module commands:

#INFILE="chr22.dose.vcf.gz"
INFILE="airwave_b37"
#OUTFILE="chr22_Airwave_CPMG_Plasma"
OUTFILE="airwave_b37"
SAMPLE_IDs="Airwave_CPMG_Plasma.txt_sample_names_FID_IID_2.txt"
SLOTS="16"

#File management:
cp $PBS_O_WORKDIR/$INFILE $TMPDIR
cp $PBS_O_WORKDIR/$SAMPLE_IDs $TMPDIR
##cp $PBS_O_WORKDIR/$CONFIG $TMPDIR

#Command:
/home/aberlang/bin/plink_1.90_3.37_16_May.dir/plink --vcf $INFILE --double-id --make-bed --keep $SAMPLE_IDs --out $OUTFILE
#/home/aberlang/bin/plink_1.90_3.37_16_May.dir/plink --bfile airwave_b37 --make-bed --keep Airwave_CPMG_Plasma.txt_sample_names_FID.txt --out airwave_b37_metabolomics_subset

#Notify by email:
##PBS -m abe


#File management:
cp $OUTFILE* $PBS_O_WORKDIR

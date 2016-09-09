#!/bin/sh

#Job parameters:
#PBS -l walltime=01:00:00
#PBS -l select=1:ncpus=16:mem=20gb
#PBS -q med-bio 

#Use current working directory
##$ -cwd

#Save standard error and out to files:
##$ -e stderr.file
##$ -o stdout.file

#Module commands:
module load gatk/3.5

INFILE="chr22.dose.vcf.gz"
OUTFILE="chr22_Airwave_CPMG_Plasma.vcf"
INFO_FILE="chr22.info.gz"
#SAMPLE_IDs="Airwave_CPMG_Plasma.txt_sample_names_FID_IID.txt"
#REFERENCE="chr22.fa"
#REFERENCE_INDEXES="chr22."
SLOTS="16"
##CONFIG=.ncbirc

#File management:
cp $PBS_O_WORKDIR/$INFILE $TMPDIR
cp $PBS_O_WORKDIR/$INFO_FILE $TMPDIR
#cp $PBS_O_WORKDIR/$REFERENCE $TMPDIR
#cp $PBS_O_WORKDIR/$REFERENCE_INDEXES* $TMPDIR
#cp $PBS_O_WORKDIR/$SAMPLE_IDs $TMPDIR
##cp $PBS_O_WORKDIR/$CONFIG $TMPDIR

#Command:
/home/aberlang/bin/gcta64.dir/gcta64 --dosage-mach-gz $INFILE $INFO_FILE --make-bed --out $OUTFILE

#Notify by email:
##PBS -m abe


#File management:
cp $OUTFILE $PBS_O_WORKDIR

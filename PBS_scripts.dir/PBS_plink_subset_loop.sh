#!/bin/sh

#Job parameters:
#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=16:mem=10gb
#PBS -q med-bio 

#Use current working directory
##$ -cwd

#Save standard error and out to files:
#PBS -e plink_loop_stderr.file
#PBS -o plink_loop_stdout.file

#Module commands:

#INFILE="chr22.dose.vcf.gz"
#OUTFILE="chr22_Airwave_CPMG_Plasma"
SAMPLE_IDs="Airwave_CPMG_Plasma.txt_sample_names_FID_IID_2.txt"
SLOTS="16"

#File management:
cp $PBS_O_WORKDIR/$SAMPLE_IDs $TMPDIR
##cp $PBS_O_WORKDIR/$CONFIG $TMPDIR

#Command:
for f in *.dose.vcf.gz; do
	cp $PBS_O_WORKDIR/$f $TMPDIR
	/home/aberlang/bin/plink_1.90_3.37_16_May.dir/plink --vcf $f \
        	                                            --double-id \
                	                                    --make-bed \
                        	                            --keep $SAMPLE_IDs \
                                	                    --out ${f}.subset
	cp *.subset* $PBS_O_WORKDIR
done

#Notify by email:
##PBS -m abe

#File management:
#cp *.subset* $PBS_O_WORKDIR

#!/bin/sh

#Job parameters:
#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=16:mem=8gb

#Use current working directory
##$ -cwd

#Save standard error and out to files:
##$ -e stderr.file
##$ -o stdout.file

#Module commands:
module load plink/1.90

INFILE=chr22.dose.vcf.gz
PHENO_FILE=screen_data_processed_plink.tsv
PHENO_NAME=C_REACTIVE_PROTEIN
COVAR_NAMES=HOURS_SINCE_EAT_BLOOD, HAS_DRUNK_WITHIN_24HOURS, IS_MENSES, IS_PREGNANT, IS_SMOKER, DIAG_STROKE, DIAG_DIABETES, DIAG_DIABETES_TYPE_1, DIAG_DIABETES_TYPE_2, DIAG_HEARTATTACK, DIAG_HIGH_CHOLESTEROL, DIAG_HYPERTENSION, DIAG_OTHER_CONDITION, BODY_MASS_INDEX, TAKES_INTENSE_EXERCISE, FAT_PERCENTAGE, age
OUTFILE=chr22_CRP
SLOTS=16
##CONFIG=.ncbirc

#File management:
cp $PBS_O_WORKDIR/$INFILE $TMPDIR
cp $PBS_O_WORKDIR/$PHENO_FILE $TMPDIR
cp $PBS_O_WORKDIR/$PHENO_NAME $TMPDIR
cp $PBS_O_WORKDIR/$COVAR_NAMES $TMPDIR
##cp $PBS_O_WORKDIR/$CONFIG $TMPDIR

#Command:
#Quantitative association analysis:

plink --noweb --vcf $INFILE \
	--pheno $PHENO_FILE \
	--pheno-name $PHENO_NAME \
	--linear sex hide-covar --adjust --ci 0.95 \
	--missing-phenotype -9 \
	--covar $PHENO_FILE \
	--covar-name $COVAR_NAMES \
	--out $OUTFILE

#File management:
cp $OUTFILE $PBS_O_WORKDIR

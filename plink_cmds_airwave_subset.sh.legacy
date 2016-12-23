# General plink commands to process Airwave data

# Plink command to subset by individuals and convert vcf to plink binary:
# in: /groupvol/med-bio/aberlang/data.dir/Airwave.dir
# /home/aberlang/bin/plink_1.90_3.37_16_May.dir/plink --vcf $INFILE --double-id --make-bed --keep $SAMPLE_IDs --out $OUTFILE

# qsub script for one file:
qsub PBS_plink_subset.sh

# Get frequency stats:
plink --bfile chr22_Airwave_CPMG_Plasma --freq --out test_chr22.test

# Convert bfile to ped:
plink --bfile chr22_Airwave_CPMG_Plasma --recode tab --out test_chr22_Airwave_CPMG_Plasma


# Run loop to convert vcf to binary and subset individuals:
# Bash loop to process with qsub (use #Ruffus!):
qsub PBS_plink_subset_loop.sh


# Create list of subsetted files:
ls -1 *.subset.bed > list_subsetted_files.txt

# Merge all resulting plink binaries into one file:
plink --make-bed --bfile [?] --merge-list list_subsetted_files.txt --out Airwave_imputed_metabolomics_subset



# Continue with recoding and converting to QTL format required:
# TO DO: This creates cuadruple IDs xxx_xxx_xxx_xxxx (use --double-id ?)
bash 00_eQTL_genotype_process.sh chr22_Airwave_CPMG_Plasma chr22_Airwave_CPMG_Plasma.matrixQTL chr22_Airwave_CPMG_Plasma.A-transpose chr22_Airwave_CPMG_Plasma.A-transpose.matrixQTL.geno

# Run update_IDs_plink_to_MatrixEQTL.R
# Run airwave_matrixeqtl_processing.R

# Run genotype QC and check non-imputed (re-run imputation?)


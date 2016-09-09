#!/bin/bash

#Use current working directory
#$ -cwd

#Run on one processor
#$ -pe dedicated 4

# select all.q queue
#$ -q all.q

#Memory option
#-l mem=10G

#Save standard error and out to files:
#$ -e stderr_eQTl_04.file
#$ -o stdout_eQTL_04.file

#Run the job - program must have full path
/bin/bash BEST-D_association_analysis_commands.sh

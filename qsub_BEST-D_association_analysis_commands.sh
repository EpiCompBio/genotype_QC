#!/bin/bash

#Use current working directory:
#$ -cwd

#Run on n processors:
#$ -pe dedicated 2 

# select all.q queue:
#$ -q all.q

#Memory option, each thread specified above will use the amount of memory specified here:
#$ -l mem_free=25G

#Save standard error and out to files:
#$ -e stderr.file
#$ -o stdout.file

#Run the job - program must have full path:
/bin/bash BEST-D_association_analysis_commands.sh 

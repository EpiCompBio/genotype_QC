#!/bin/bash

#Use current working directory
#$ -cwd

#Run on one processor
#$ -pe dedicated 3 

# select all.q queue
#$ -q all.q

#Memory option, each thread specified above will use the amount of memory specified here:
#$ -l mem_free=25G

#Save standard error and out to files:
#$ -e plink_IBS_stderr.file
#$ -o plink_IBS_stdout.file

#Prepare a second, executable script that is launched from this one:
#Run the job - program must have full path:
./qsub_plink_IBS.sh

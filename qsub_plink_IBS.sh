#!/bin/bash

#Run the job - program must have full path:
/ifs/apps/bio/plink-1.07/plink-1.07-x86_64/plink --noweb --bfile P140343-Results_FinalReport --extract plink.prune.in --genome --out P140343-Results_FinalReport_pruned

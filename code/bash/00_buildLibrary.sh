#!/bin/bash
# Ruth Gómez Graciani
# 30 06 2021

###############################################################################
# Description:                                                                
# Build a library                 
###############################################################################

SCRIPTPATH=$1
FASTAFILE=$2


cd $SCRIPTPATH


 cp $FASTAFILE data/use/bowtie_index/
 cd data/use/bowtie_index/

 bowtie2-build bplib.fa bplib
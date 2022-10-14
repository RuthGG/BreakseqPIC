#!/bin/bash
# Ruth Gómez Graciani
# 30 06 2021

###############################################################################
# Description:                                                                
# Build a library                 
###############################################################################

SCRIPTPATH=$1
DATAPATH=$2
FASTAFILE="${DATAPATH}/bplib.fa"
COORDSFILE="${DATAPATH}/bplib.coords"

# Build library in external storage
cd $SCRIPTPATH

 cp $FASTAFILE data/use/bowtie_index/
 cd data/use/bowtie_index/

 bowtie2-build bplib.fa bplib

# Check library within analysis directory
mkdir -p "${SCRIPTPATH}/${DATAPATH}/coordCheck/"
cd  "${SCRIPTPATH}/${DATAPATH}/"

bowtie2 -f -x "${SCRIPTPATH}/data/use/bowtie_index/h_sapiens_asm" -U bplib.fa > "coordCheck/result"
Rscript ${SCRIPTPATH}/code/rscript/extra_checkCoordinates.R


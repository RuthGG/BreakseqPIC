#!/bin/bash
# Ruth Gómez Graciani
# 29 09 2022

SRACODE=$1 # code
SRADIR=$2  # absolute path from SINGULARITY HOME, includes DATE_NAME
FASTQDIR=$3 # absolute path from SINGULARITY HOME; "" to override fasterq-dump, includes DATE_NAME
PASSWORD=$4 # abolute path from SINGULARITY HOME

# Create dirs if they do not exist
mkdir -p $SRADIR $FASTQDIR

# Save current path
MYBASE=$(pwd)

cd $SRADIR

#---------------------------------- 
# Prefetch code
# ---------------------------------
if [ ! -f $SRACODE/*.sra ]; then #if SRA does not exist

 exec 3>&1 4>&2
  trap 'exec 2>&4 1>&3' 0 1 2 3
exec 1> ${SRACODE}.log 2>&1

prefetch --max-size 100000000 --ngc "${MYBASE}/${PASSWORD}" $SRACODE

fi

#---------------------------------- 
# fastq conversion
# ---------------------------------

if [[ $FASTQDIR != "" ]]; then
    echo "Converting SRA file ${SRACODE} to FASTQ"
    mkdir -p ${MYBASE}/${FASTQDIR}/${SRACODE}
    fasterq-dump --ngc "${MYBASE}/${PASSWORD}" -O ${MYBASE}/${FASTQDIR}/${SRACODE}/ ./${SRACODE}
fi

echo "Done!"


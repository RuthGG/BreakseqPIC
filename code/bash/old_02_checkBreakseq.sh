#!/bin/bash
# Ruth Gómez Graciani
# 04 06 2021

###############################################################################
# Description:                                                                
# check if breakseq worked                        
###############################################################################

DIR=$1
INDFILE=$2

SCRIPTPATH=$(pwd)

	# Check if important file is 0
	cd ${DIR}

	numbers=$(ls -la  */*.uni.sam | awk '$5 > 0{print $9}'|  sed 's/\/.*//g' | grep -f - -v -n ${SCRIPTPATH}/${INDFILE} | sed 's/:.*//' | awk '{print ($1-1)}' )
	length=$(cat  ${SCRIPTPATH}/${INDFILE}   | wc -l | wc -c | awk '{print $1-1}')

	numbers=$(for n in $numbers; do  length_num=$(echo $n | wc -c | awk '{print $1-1}' ); while [ $length_num -lt $length  ]; do n=$(echo 0$n); length_num=$(echo $n | wc -c | awk '{print $1-1}' ); done; echo $n; done)

	echo $numbers | tr ' ' ','
	echo $numbers | tr ' ' ',' > failedRefs.txt
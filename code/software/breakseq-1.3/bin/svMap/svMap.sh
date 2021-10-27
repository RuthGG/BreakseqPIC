#!/bin/sh

if [ $# -lt 3 ]
then
	echo "Usage: $0 <reference> <library> <output dirname/basename> <min_cover> <lib_len> <query>..."
	exit 0
fi

echo "## Breakpoint Library Mapping started at `date` ##"

export SEED_LEN=36 # 36 
export MIN_COVER=$4 # cambia, 15, 20...
export LIB_LEN=$5
export SCORE_MIN=$6

BIN=`cd \`dirname $0\`; pwd`
REF=$1
LIB=$2
PREFIX=$3_

shift 6

for QUERY in $*
do
	echo "## Processing query $QUERY ##"
	OUTPUT=$PREFIX`basename $QUERY`

	echo "## Aligning the query against the library ##"
	sh $BIN/align $LIB $QUERY $SCORE_MIN > $OUTPUT.ini.sam

	echo "## Filtering hits for breakpoint coverage ##"
	# [2012-10-26] Substituted by Ignasi.
	#/usr/bin/python2 $BIN/filter.py $LIB_LEN $SEED_LEN $MIN_COVER $OUTPUT.ini.sam > $OUTPUT.fil.sam
	perl $BIN/filter.plx $LIB_LEN $SEED_LEN $MIN_COVER $OUTPUT.ini.sam > $OUTPUT.fil.sam

	echo "## Converting hits to Fastq ##"
	/usr/bin/python2 $BIN/sam2fq.py $OUTPUT.fil.sam > $OUTPUT.fil.fq

	echo "## Aligning the hits against the reference ##"
	sh     $BIN/align $REF $OUTPUT.fil.fq $SCORE_MIN > $OUTPUT.ref.sam 

	echo "## Splitting hits for uniqueness ##"
	/usr/bin/python2 $BIN/splituni.py $OUTPUT.fil.sam $OUTPUT.ref.sam > $OUTPUT.uni.sam 2> $OUTPUT.xun.sam
done


echo "## Scoring and reporting calls ##"
SVCALL=$PREFIX`date +%y%m%d_%H%M%S`.svcall
/usr/bin/python2 $BIN/svhit.py $PREFIX*.uni.sam > $SVCALL.uni
/usr/bin/python2 $BIN/svhit.py $PREFIX*.xun.sam > $SVCALL.xun
/usr/bin/python2 $BIN/svcall.py $SVCALL.uni $SVCALL.xun > $SVCALL 2> $SVCALL.txt

echo "## Breakpoint Library Mapping ended at `date` ##"

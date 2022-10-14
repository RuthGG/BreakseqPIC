#!/bin/bash
#
# This needs the bam files.
#
# It will produce the following files:
#
#   - genotypes (not overwritten)
#   - filtered.uni (overwritten)
#   - filtered.xun (overwritten)
#   - filtered.results (overwritten)
#   - filtered.results.genotyped (overwritten)
#   - filtered.errors (overwritten)
#
# After processing of the mapped reads, the available evidence
# to predict the genotypes is filtered, cleaner, and it is now that the genotypes
# should be predicted, using uni.sorted.bam and xun.sorted.bam, if available.
#

#if [ $# -lt 1 ]; then
#	echo "Usage FilterUniXun.sh <dir>"
#	echo "     <dir>: Current directory, where individual's folders must be."
#	exit
#fi


#########
## CPU ##
#########
CPU=1
FOLDER_PRUEBA=$1

##########################################
## CAMBIAR DIRECCIONES SEGUN CONVENGA!! ##
##########################################
#DATADIR=/home/jon/Desktop/Functional_Analysis/HsInv0102/CLL_genotyping/breakseq_analysis/datos_hsinv0102
#DATADIR=/home/jon/Desktop/Functional_Analysis/HsInv0102/Breakseq_1KGP/Results_allINVs/datos
#DATADIR=/home/jon/Desktop/Functional_Analysis/HsInv0102/Breakseq_1KGP/Results_allINVs/datos_newName
#PYTHON=/usr/local/bin/python
# PYTHON=python2
SVHIT=/breakseq-1.3/bin/svMap/svhit.py
BIN=../../code/aux_breakseq
#WORKDIR=/home/jon/Desktop/Functional_Analysis/HsInv0102/CLL_genotyping/breakseq_analysis/results
#WORKDIR=/home/jon/Desktop/Functional_Analysis/HsInv0102/CLL_metaanalysis/results/HsInv0102
#WORKDIR=/home/jon/Desktop/Functional_Analysis/HsInv0102/Breakseq_1KGP/Results_allINVs
#WORKDIR=/home/jon/Desktop/Functional_Analysis/HsInv0102/Breakseq_1KGP/Results

## Analisis 102 - metaanalysis2 - rutas
# DATADIR=/home/jon/Desktop/Functional_Analysis/HsInv0102/CLL_genotyping/breakseq_analysis/datos_hsinv0102
# WORKDIR=/home/jon/Desktop/Functional_Analysis/HsInv0102/CLL_metaanalysis2/results/HsInv0102

# ## St Jude - rutas
# DATADIR=/home/jon/Desktop/Functional_Analysis/HsInv0102/CLL_genotyping/breakseq_analysis/datos_hsinv0102
# WORKDIR=/home/jon/Desktop/Functional_Analysis/HsInv0102/StJude/Step4_Breakseq/HsInv0102/

## Analysis allINVs (1KGP) - Paper funcional
DATADIR=../../data/raw/datos_librerias
WORKDIR=./



for i in `cat $WORKDIR/individuals.txt`; do
	if [ -e $i/uni.sorted.bam ] && [ ! -e $i/filtered.uni ]; then
		samtools view $i/uni.sorted.bam > $i/uni.sorted.sam
		# head -n-3 $i/uni.sorted.sam > $i/uni.sorted.sam
		python2 /breakseq-1.3/bin/svMap/svhit.py $i/uni.sorted.sam > $i/filtered.uni
		# rm $i/uni.sorted.sam
	fi
	if [ -e $i/xun.sorted.bam ] && [ ! -e $i/filtered.xun ]; then
		samtools view $i/xun.sorted.bam > $i/xun.sorted.sam
		python2 /breakseq-1.3/bin/svMap/svhit.py $i/xun.sorted.sam > $i/filtered.xun
		# rm $i/xun.sorted.sam
	fi
done

#for i in `find . -name 'filtered.xun'`; do
#	S=`awk -v I=$i 'BEGIN{split(I,A,/\//); print A[2]}'`
#       I=`awk -v I=$i 'BEGIN{sub(/.xun$/,"",I); print I}'`
#       ../../bin/result_table.2.pl -i $I -s $S >> z1
#done



# ORIGINAL
# for i in `find . -type d`; do
#  	S=`awk -v I=$i 'BEGIN{split(I,A,/\//); print A[2]}'`
#  	$BIN/result_table.2_allINVs_cpuN.pl -i $i/filtered -s $S >> z1
# done

# NEW
for i in `ls -d */`; do
	# NEW
	S=`echo $i | sed 's/\///'`
	# $BIN/result_table.2.pl -i $S/filtered -s $S >> z1
	#  ORIGINAL
	# S=`awk -v I=$i 'BEGIN{split(I,A,/\//); print A[2]}'`
 	$BIN/result_table.2_allINVs_cpuN.pl -i $i/filtered -s $S >> z1
done

sort -k 1g,1 -k 4,4 -k 2,2 -k 3,3 z1 > filtered.results

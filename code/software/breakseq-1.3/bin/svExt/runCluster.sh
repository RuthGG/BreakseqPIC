#!/bin/sh

if [ $# -lt 1 ]
then
	echo "Please specify the SVs input in a pre-defined txt format"
	exit 0
fi

work=`cd \`dirname $0\`; pwd`
base=`cd $work/../../; pwd`

input=$1
prefix=${input/.txt/}
gff=$prefix.gff

echo "** Converting txt into gff **"
$work/txt2gff.py $input > $gff

echo "** Clustering input by different size groups  **"
$work/svSize.py $gff `dirname $input` 50 200
$work/svSize.py $gff `dirname $input` 200 1000
$work/svSize.py $gff `dirname $input` 1000

echo "** Submitting jobs"
timestamp=`date +%Y%m%d-%H%M`
qin "$base/breakseq annotate $prefix.50-200bp.gff   $timestamp/sv50-200bp   50"
qin "$base/breakseq annotate $prefix.200-1000bp.gff $timestamp/sv200-1000bp 200"
qin "$base/breakseq annotate $prefix.1000-maxbp.gff $timestamp/sv1000-maxbp 1000"

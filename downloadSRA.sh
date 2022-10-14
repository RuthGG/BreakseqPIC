#!/bin/bash
# Ruth Gómez Graciani
# 29 09 2022

DATE=$(date +%F)

SRADIR=$1

SRALIST=${SRADIR}/accession_download.txt
SRAINVS=${SRADIR}/indiv_list.txt
SRAEQ=${SRADIR}/accession_list.txt
PASSWORD=`ls ${SRADIR}/*.ngc` # from inside 20210325_breakseq

waitForCondor () {
	while [[ $(condor_q | condor_q  | grep Total |  head -n 2 | tail -n 1| awk '{print $10}') -gt 0 ]] || [[ $(condor_q | condor_q  | grep Total |  head -n 2 | tail -n 1 | awk '{print $12}') -gt 0 ]]; do
	
		condor_q
		sleep ${1}m

	done

	echo "All done with Condor! Jobs on hold:"
	condor_q | condor_q  | grep Total |  head -n 2 | tail -n 1 | awk '{print $14}'

}

mkdir -p /data/bioinfo/common/breakseq_data/${DATE}_SRA_avery

# SRACODE=$1 # code
# SRADIR=$2  # absolute path from SINGULARITY HOME
# FASTQDIR=$3 # absolute path from SINGULARITY HOME; "" to override fasterq-dump
# PASSWORD=$4 # abolute path from SINGULARITY HOME

	echo "executable = /bin/singularity
args = \"exec \\
--bind /data/bioinfo/common/breakseq_data:/nfs/pic.es/user/r/rgomez/20210325_breakseq/tmp/downloadSRA:rw \\
/data/bioinfo/software/sratoolkit.2.11.3_latest.sif  \\
bash 20210325_breakseq/code/bash/01_downloadSRA.sh \$(item) '20210325_breakseq/tmp/downloadSRA/${DATE}_SRA_avery' '' '20210325_breakseq/$PASSWORD' \"

output = tmp/condor.out
error = tmp/condor.err
log = tmp/condor.log

request_cpus = 1

queue 1 from ${SRALIST}" > tmp/downloadSRA.sub

condor_submit tmp/downloadSRA.sub

waitForCondor 15

# mkdir /data/bioinfo/scratch/breakseq_fastqs/${DATE}_fastq_avery/
# for i in `cat $SRAINVS` ; do
# 	echo $i
# 	LIST=`grep "${i}[a-z]*\s" data/raw/avery_data/accession_list.txt | cut -f2`
# 	LISTGREP=`echo $LIST | sed 's/ /\|/g' `
# 	echo $LISTGREP 
# 	FILELIST=`ls /data/bioinfo/scratch/breakseq_fastqs/${DATE}_SRA_avery/*/*.fastq | grep -E "$LISTGREP"`
# 	echo $FILELIST
# 	mkdir -p /data/bioinfo/scratch/breakseq_fastqs/${DATE}_fastq_avery/$i
# 	if [[ $FILELIST != "" ]]; then
# 		cat $FILELIST > /data/bioinfo/scratch/breakseq_fastqs/${DATE}_fastq_avery/$i/selected_regions.fastq
# 		rm $FILELIST
# 	fi
# 	FILELIST=""
# 	LIST=""
# done

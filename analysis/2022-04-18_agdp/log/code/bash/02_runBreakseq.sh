#!/bin/bash
# Ruth Gómez Graciani
# 14 05 2021

###############################################################################
# Description:                                                                
# Run breakseq analysis                            
###############################################################################

# GET VARIABLES - there are examples commented
# =========================================================================== #
TMPDIR=$1     
OUTDIR=$2  
FASTQ_DIR=$3   
MIN_COVER=$4
LIB_LEN=$5
SCORE_MIN=$6
SAMPLE=$7

# SET DERIVED VARIABLES 
# =========================================================================== #

#  RUN BREAKSEQ
# =========================================================================== #

echo "# Running BreakSeq on $SAMPLE."

# Make a folders for the sample
mkdir -p ${OUTDIR}/${SAMPLE} ${TMPDIR}/${SAMPLE}
    
# Run breakseq
./code/software/breakseq-1.3/breakseq identify ${TMPDIR}/${SAMPLE}/ $MIN_COVER $LIB_LEN $SCORE_MIN ${FASTQ_DIR}/${SAMPLE}/selected_regions.fastq  1> ${TMPDIR}/${SAMPLE}/bpseq.log 2> ${TMPDIR}/${SAMPLE}/bpseq.err
        
echo "Breakseq log"
cat ${TMPDIR}/${SAMPLE}/bpseq.log
echo "Breakseq err"
cat ${TMPDIR}/${SAMPLE}/bpseq.err

mv ${TMPDIR}/${SAMPLE}/*ref.sam ${TMPDIR}/${SAMPLE}/*uni.sam ${TMPDIR}/${SAMPLE}/*xun.sam ${OUTDIR}/${SAMPLE}
cat ${TMPDIR}/${SAMPLE}/*ini.sam  | awk '{a[$3]+=1} END{for (i in a) print i,a[i]}' >> ${OUTDIR}/${SAMPLE}/inisam_summary
cat ${TMPDIR}/${SAMPLE}/*fil.sam  | awk '{a[$3]+=1} END{for (i in a) print i,a[i]}' >> ${OUTDIR}/${SAMPLE}/filsam_summary

# Add rm tmp
rm -r ${TMPDIR}/${SAMPLE}/
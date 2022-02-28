#!/bin/bash
# Ruth Gómez Graciani
# 29 06 2021

###############################################################################
# Description:                                                                
# Top-level run of commands; creates recipes, parallelizes and summarizes        
###############################################################################

## Recieve parameters ## 
# =========================================================================== #

DATE=$(date +%F)
source runCommands_config.sh

## Set ennvironment ##
# =========================================================================== #

OUTDIR="analysis/${DATE}_${NAME}"
TMPDIR="tmp/${DATE}_${NAME}"
DATADIR="analysis/${DATE}_${NAME}/data"
LOGDIR="analysis/${DATE}_${NAME}/log"
BASEPATH=$(pwd)

mkdir -p $TMPDIR $DATADIR $LOGDIR

cp runCommands_config.sh ${LOGDIR}/runCommands_config.sh

## Set functions ##
# =========================================================================== #

waitForCondor () {
	while [[ $(condor_q | condor_q  | grep Total |  head -n 2 | tail -n 1| awk '{print $10}') -gt 0 ]] || [[ $(condor_q | condor_q  | grep Total |  head -n 2 | tail -n 1 | awk '{print $12}') -gt 0 ]]; do
	
		condor_q
		sleep ${1}m

	done

	echo "All done with Condor! Jobs on hold:"
	condor_q | condor_q  | grep Total |  head -n 2 | tail -n 1 | awk '{print $14}'

}

checkDownload () {

		# How many invs do we want to download (plus unmapped)?
		TARGETS=$(( $(cat $INVSFILE | wc -l)+1 ))

		# Make a list of fails
		# This is the list for the queue
		wc -l ${TMPDIR}/01_download/readscount/* | sed 's/^ *//g' | awk  -v t="$TARGETS" '$1==t{print $0}' | cut -d " " -f2 | sed 's/^.*readscount\///g' | sed 's/\.txt$//g'  > ${TMPDIR}/01_download/successnames

		grep -v -f ${TMPDIR}/01_download/successnames $SAMPLESLIST > ${TMPDIR}/01_download/failednames
		# Counter
		TOTEST=$(cat ${TMPDIR}/01_download/failednames | wc -l)

		echo "Downloads tested. $TOTEST inds with errors, about to be repeated"

}

echo "# STEP 00 - Library
# =========================================================================== #"


if [[ $OVERRIDE -gt 0 ]]; then
	echo "This step will not run now"

else

	mkdir -p ${LOGDIR}/00_library/ ${DATADIR}/datos_librerias 

	cp data/raw/seed_librerias/${LIBRARY_SEED}/bplib* data/raw/seed_librerias/${LIBRARY_SEED}/ref* ${DATADIR}/datos_librerias
	cd ${DATADIR}/datos_librerias

	# Generate probe header
	# CHANGE IF WE CHANGE LIBRARY LENGTH!!
	grep '>' bplib.fa | sed -e 's/>//g' | sed -e 's/^/@SQ\tSN:/g' | sed -e "s/$/\tLN:$LIBRARY_LENGTH/g" > bplib.header.template.txt

	echo '@HD VN:1.0 SO:unsorted' | tr ' ' '\t' > tmp
	echo '@RG ID:XXX SM:XXX' | tr ' ' '\t' >> tmp
	echo '@PG ID:breakseq' | tr ' ' '\t' >> tmp
	cat bplib.header.template.txt >> tmp
	rm bplib.header.template.txt
	mv tmp bplib.header.template.txt

	# Generate inversion list 
	cat  bplib.coords | cut -f1 | sed 's/^\w\w\w//'|sed 's/BP.*$//' | sort | uniq > inversions_completeList.txt

	# Make recipe to change current index

	cd ${BASEPATH}/

	echo "executable = /bin/singularity
args = \"exec \\
--bind /data/bioinfo/common/bowtie2_index:/nfs/pic.es/user/r/rgomez/20210325_breakseq/data/use/bowtie_index:rw \\
/data/bioinfo/software/rgomez_breakseq.sif \\
bash 20210325_breakseq/code/bash/00_buildLibrary.sh ${BASEPATH} ${DATADIR}/datos_librerias/\"

output = ${LOGDIR}/00_library/condor.out
error = ${LOGDIR}/00_library/condor.err
log = ${LOGDIR}/00_library/condor.log

request_cpus = 1

queue 1" > ${LOGDIR}/00_library.sub

	if [[ $LIBRARY_UPDATE == "yes" ]]; then

		condor_submit ${LOGDIR}/00_library.sub

		waitForCondor 1 

		cp -r ${DATADIR}/datos_librerias/coordCheck data/raw/seed_librerias/${LIBRARY_SEED}/
	fi
fi


## This will always run, it sets parameters for later ##

if [[ -z $INVSFILE ]]; then
	INVSFILE="${DATADIR}/datos_librerias/inversions_completeList.txt"
fi

	# This is to avoid touching raw data files + keeping log about the process
	cp $INVSFILE "${DATADIR}/regions.txt"
	cat $SAMPLESLIST | cut -f1 > "${DATADIR}/samples.txt"
	cp $SAMPLESFILE "${DATADIR}/samples_pathIndex.txt"

	INVSFILE="${DATADIR}/regions.txt"
	SAMPLESLIST="${DATADIR}/samples.txt"
	SAMPLESFILE="${DATADIR}/samples_pathIndex.txt"


#######################################################


echo "# Step 01 - Download
# =========================================================================== #"

if [[ $OVERRIDE -gt 1 ]] || [[ ! -z $READ_FASTQS ]]; then
	echo "This step will not run now"
	mkdir -p ${OUTDIR}/01_download/
	cp ${READ_FASTQS}_readscount.txt ${OUTDIR}/01_download/readscount.txt
else

	mkdir -p ${LOGDIR}/01_download/ ${OUTDIR}/01_download/ ${TMPDIR}/01_download/ ${DATADIR}/bamFiles/
	
	# make tmp file
	if [[ $RESUME == "y" ]]; then
		echo "RESUMING donwload"
		checkDownload 
	else
		echo "STARTING download (data replacement may occur)"
		cp $SAMPLESLIST ${TMPDIR}/01_download/failednames
		TOTEST=1
	fi

	# make recipe
	echo "executable = /bin/singularity
args = \"exec \\
--bind /data/bioinfo/scratch/breakseq_fastqs:/nfs/pic.es/user/r/rgomez/20210325_breakseq/${OUTDIR}/01_download/:rw \\
--bind /data/bioinfo/scratch/breakseq_bam:/nfs/pic.es/user/r/rgomez/20210325_breakseq/${DATADIR}/bamFiles/:rw \\
/data/bioinfo/software/rgomez_breakseq.sif \\
bash 20210325_breakseq/runall.sh download -r ${INVSFILE} -s \$(item) -t \$(item) -n ${DATE}_${NAME}\"

output = ${LOGDIR}/01_download/condor.out
error = ${LOGDIR}/01_download/condor.err
log = ${LOGDIR}/01_download/condor.log

request_cpus = 1

queue 1 from ${TMPDIR}/01_download/failednames" > ${LOGDIR}/01_download.sub

	# set loop
	
	LOOPCOUNT=0
	while [[ $TOTEST -gt 0 ]]; do
		LOOPCOUNT=$(($LOOPCOUNT+1))
		echo "Starting loop $LOOPCOUNT"

		condor_submit ${LOGDIR}/01_download.sub
	
		waitForCondor 15 

		# CHECK CORRECT DOWNLOAD
		checkDownload 
	done

	rm ${TMPDIR}/01_download/failednames ${TMPDIR}/01_download/successnames

	# Make summary counts and copy them in folders
	cat ${TMPDIR}/01_download/readscount/* > ${OUTDIR}/01_download/readscount.txt
	cp ${OUTDIR}/01_download/readscount.txt /data/bioinfo/scratch/breakseq_fastqs/${DATE}_${NAME}_readscount.txt
	

fi

## This will always run, it sets parameters for later ##

if [[ -z $READ_FASTQS ]]; then
	READ_FASTQS=/data/bioinfo/scratch/breakseq_fastqs/${DATE}_${NAME}/
	mkdir -p ${OUTDIR}/01_download/
	cp /data/bioinfo/scratch/breakseq_fastqs/${DATE}_${NAME}_readscount.txt ${OUTDIR}/01_download/readscount.txt
fi

#######################################################

echo "# Step 02 - Breakseq
# =========================================================================== #"

if [[ $OVERRIDE -gt 2 ]] ; then
	echo "This step will not run now"
else

	mkdir -p ${LOGDIR}/02_breakseq/ ${OUTDIR}/02_breakseq/ ${TMPDIR}/02_breakseq
	
	# make tmp file
	cp $SAMPLESLIST ${TMPDIR}/02_breakseq/failednames


	echo "executable = /bin/singularity
args = \"exec \\
--bind ${READ_FASTQS}:/nfs/pic.es/user/r/rgomez/20210325_breakseq/${OUTDIR}/01_download/:rw \\
--bind /data/bioinfo/common/bowtie2_index:/nfs/pic.es/user/r/rgomez/20210325_breakseq/data/use/bowtie_index \\
--bind /data/bioinfo/scratch/breakseq_tmp/:/nfs/pic.es/user/r/rgomez/20210325_breakseq/${TMPDIR}/:rw \\
/data/bioinfo/software/rgomez_breakseq.sif \\
bash 20210325_breakseq/runall.sh breakseq -k $SCORE_MIN -m $COV_AROUND -l $LIBRARY_LENGTH -f ${OUTDIR}/01_download/ -s \$(Item) -t \$(Item) -n ${DATE}_${NAME}\"

output = ${LOGDIR}/02_breakseq/condor.out
error = ${LOGDIR}/02_breakseq/condor.err
log = ${LOGDIR}/02_breakseq/condor.log

request_cpus = 1

queue 1 from ${TMPDIR}/02_breakseq/failednames "  >  ${LOGDIR}/02_breakseq.sub
	
	# set loop
	TOTEST=1
	LOOPCOUNT=0
	while [[ $TOTEST -gt 0 ]]; do
		LOOPCOUNT=$(($LOOPCOUNT+1))
		echo "Starting loop $LOOPCOUNT"


		condor_submit ${LOGDIR}/02_breakseq.sub
	
		waitForCondor 15 

		# CHECK CORRECT DOWNLOAD
		# -----------------------------------------------------

		# how many inds were downloaded?
		ls -l ${OUTDIR}/02_breakseq/*/*uni.sam | awk '$5>0{print $9}' | sed "s/^.*02_breakseq\///g" | sed 's/\/.*//g' > ${TMPDIR}/02_breakseq/successnames

		# Make a list of fails
		# This is the list for the queue
		grep -v -f ${TMPDIR}/02_breakseq/successnames $SAMPLESLIST  > ${TMPDIR}/02_breakseq/failednames

		# Counter
		TOTEST=$(cat ${TMPDIR}/02_breakseq/failednames | wc -l)


	done

	rm ${TMPDIR}/02_breakseq/failednames ${TMPDIR}/02_breakseq/successnames
	cat  ${OUTDIR}/02_breakseq/*/inisam_summary > ${OUTDIR}/02_breakseq/inisam_summary
	cat  ${OUTDIR}/02_breakseq/*/filsam_summary > ${OUTDIR}/02_breakseq/filsam_summary

	rm ${OUTDIR}/02_breakseq/*/inisam_summary ${OUTDIR}/02_breakseq/*/filsam_summary




fi


## This will always run, it sets parameters for later ##

if [[ -z $BREAKSEQ_RESULTS ]]; then
	BREAKSEQ_RESULTS="$OUTDIR/02_breakseq/"
fi

#######################################################



echo "# Step 03 - ProcessAligned
# =========================================================================== #"

if [[ $OVERRIDE -gt 3 ]] ; then
	echo "This step will not run now"
else

	mkdir -p ${LOGDIR}/03_processaligned/ ${OUTDIR}/03_processaligned/

	echo "#!/bin/bash
INREF='${INREF}'
NOTINREF='${NOTINREF}'" > $DATADIR/03_processaligned_parameters.sh
	
	echo "executable = /bin/singularity
args = \"exec \\
--bind /data/bioinfo/common/bowtie2_index:/nfs/pic.es/user/r/rgomez/20210325_breakseq/data/use/bowtie_index \\
--bind /data/bioinfo/scratch/breakseq_tmp/:/nfs/pic.es/user/r/rgomez/20210325_breakseq/${TMPDIR}/:rw \\
/data/bioinfo/software/rgomez_breakseq.sif \\
bash 20210325_breakseq/runall.sh processaligned -a $BREAKSEQ_RESULTS -s $SAMPLESLIST -m \$(Item) -t min_${COV_AROUND}_\$(Item) -n ${DATE}_${NAME}\"

output = ${LOGDIR}/03_processaligned/condor.out
error = ${LOGDIR}/03_processaligned/condor.err
log = ${LOGDIR}/03_processaligned/condor.log

request_cpus = 1

queue 1 in ($MIN_LENGTH) " >  ${LOGDIR}/03_processaligned.sub

	condor_submit ${LOGDIR}/03_processaligned.sub
		
	waitForCondor 5 

fi


## This will always run, it sets parameters for later ##

if [[ -z $GENOTYPES_FILE ]]; then
	GENOTYPES_DIR="$OUTDIR/03_processaligned/"
fi

#######################################################


echo "# Step 04 - QualityAnalysis
# =========================================================================== #"

if [[ $OVERRIDE -gt 4 ]] ; then
	echo "This step will not run now"
else

	mkdir -p ${LOGDIR}/04_qualityanalysis/ ${OUTDIR}/04_qualityanalysis/ ${TMPDIR}/04_qualityanalysis/
	
	echo "executable = /bin/singularity
args = \"exec \\
--bind /data/bioinfo/scratch/1KGP_hg19_raw/:/nfs/pic.es/user/r/rgomez/20210325_breakseq/data/raw/1KGP_data/vcf \\
--bind /data/bioinfo/scratch/breakseq_tmp/:/nfs/pic.es/user/r/rgomez/20210325_breakseq/${TMPDIR}/:rw \\
/data/bioinfo/software/rgomez_plots.sif \\
bash 20210325_breakseq/runall.sh qualityanalysis \\
-g $GENOTYPES_DIR -s $SAMPLESLIST -r ${INVSFILE} -m $MAX_ERROR -t min_${COV_AROUND}_${MIN_LENGTH} -n ${DATE}_${NAME}\"

output = ${LOGDIR}/04_qualityanalysis/condor.out
error = ${LOGDIR}/04_qualityanalysis/condor.err
log = ${LOGDIR}/04_qualityanalysis/condor.log

request_cpus = 1

queue 1" >  ${LOGDIR}/04_qualityanalysis.sub

	condor_submit ${LOGDIR}/04_qualityanalysis.sub
		
	waitForCondor 15 

fi






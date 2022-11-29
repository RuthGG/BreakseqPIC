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
export APPTAINER_LOCALCACHEDIR=/data/bioinfo/scratch/breakseq_tmp/singularityCache/
export APPTAINER_TMPDIR=/data/bioinfo/scratch/breakseq_tmp/singularityTmp/

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

	cp data/raw/seed_librerias/${LIBRARY_SEED}/bplib.* data/raw/seed_librerias/${LIBRARY_SEED}/ref* data/raw/seed_librerias/${LIBRARY_SEED}/specs.sh ${DATADIR}/datos_librerias 
	source data/raw/seed_librerias/${LIBRARY_SEED}/specs.sh

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
--bind /data/bioinfo/common/${ASSEMBLY}:/nfs/pic.es/user/r/rgomez/20210325_breakseq/data/use/bowtie_index:rw \\
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
	echo "This step will not run now. Downloads will not be deleted"
	mkdir -p ${OUTDIR}/01_download/
	# cp ${READ_FASTQS}_readscount.txt ${OUTDIR}/01_download/readscount.txt
	KEEP_DOWNLOADS="y"
else
	> /data/bioinfo/scratch/breakseq_fastqs/${DATE}_${NAME}_readscount.txt
	mkdir -p ${LOGDIR}/01_download/ ${OUTDIR}/01_download/ ${TMPDIR}/01_download/ ${DATADIR}/bamFiles/
	
	DONWSELECTOR=$(head -n1 $SAMPLESFILE | rev| cut -d "." -f1 |rev)
	
	if [[ $DONWSELECTOR == "sra" ]]; then
		
		# Take SRA list from pathindex
		grep -f ${SAMPLESLIST} ${SAMPLESFILE} | rev | cut -d"/" -f2 | rev > ${DATADIR}/sra_download.txt 
		SRALIST=${DATADIR}/sra_download.txt 
		ANCHOR=$(head -n1 ${SRALIST})

		# Take paths
		SRAPATH=$(head -n1 pathIndex.txt | cut -f2 | sed "s/${ANCHOR}.*//g") # this one will point to bamfiles
		FASTQPATH=/data/bioinfo/scratch/breakseq_tmp/${DATE}_${NAME}/ # this one will point to 01_download
		mkdir -p $FASTQPATH

		# Trasnform those SRAs into fastq (downloaded as needed)
		# SRACODE=$1 # code
		# SRADIR=$2  # absolute path from SINGULARITY HOME, includes DATE_NAME
		# FASTQDIR=$3 # absolute path from SINGULARITY HOME, , includes DATE_NAME; "" to override fasterq-dump
		# PASSWORD=$4 # abolute path from SINGULARITY HOME

		echo "executable = /bin/singularity
args = \"exec \\
--bind ${SRAPATH}:/nfs/pic.es/user/r/rgomez/20210325_breakseq/${DATADIR}/bamFiles/:rw \\
--bind ${FASTQPATH}:/nfs/pic.es/user/r/rgomez/20210325_breakseq/${OUTDIR}/01_download/:rw \\
/data/bioinfo/software/sratoolkit.2.11.3_latest.sif  \\
bash 20210325_breakseq/code/bash/01_downloadSRA.sh \$(item) 20210325_breakseq/${DATADIR}/bamFiles/ 20210325_breakseq/${OUTDIR}/01_download/ 20210325_breakseq/$NGC_PATH \"

output = ${LOGDIR}/01_download/condor.out
error = ${LOGDIR}/01_download/condor.err
log = ${LOGDIR}/01_download/condor.log

request_cpus = 1

queue 1 from ${SRALIST}" > ${LOGDIR}/01_download.sub
		
		condor_submit ${LOGDIR}/01_download.sub

		waitForCondor 15
	
		# Aggregate fastqs by the individual
		mkdir /data/bioinfo/scratch/breakseq_fastqs/${DATE}_${NAME}
		> ${LOGDIR}/01_download/fastq_aggregate
		for i in `cat $SAMPLESLIST` ; do
			echo $i >> ${LOGDIR}/01_download/fastq_aggregate
			LIST=`grep $i ${SAMPLESFILE} | rev | cut -d"." -f2 |cut -d"/" -f1 |rev`
			LISTGREP=`echo $LIST | sed 's/ /\|/g' `
			echo $LISTGREP  >> ${LOGDIR}/01_download/fastq_aggregate
			FILELIST=`ls FASTQPATH/*/*.fastq | grep -E "$LISTGREP"`
			echo $FILELIST  >> ${LOGDIR}/01_download/fastq_aggregate
			mkdir -p /data/bioinfo/scratch/breakseq_fastqs/${DATE}_${NAME}/$i
			if [[ $FILELIST != "" ]]; then
				>/data/bioinfo/scratch/breakseq_fastqs/${DATE}_${NAME}/$i/selected_regions.fastq
				for THISFILE in $FILELIST; do
					cat $THISFILE >> /data/bioinfo/scratch/breakseq_fastqs/${DATE}_${NAME}/$i/selected_regions.fastq
					rm $THISFILE
				done
			fi
			FILELIST=""
			LIST=""
		done 

	else #THIS IS FOR NORMAL AND 1KGP

		# make tmp file
		if [[ $RESUME == "y" ]]; then
			echo "RESUMING download (absent/incomplete samples from 0)" # A failed regions should exist but it will be incomplete because maybe we interrupted the download; resumed individuals will resume from start
			checkDownload 
		else
			echo "STARTING download (all data from 0)"
			# Make failednames
			cp $SAMPLESLIST ${TMPDIR}/01_download/failednames
			TOTEST=1
		fi

		# Make a base failedregions list
		> ${TMPDIR}/01_download/failedregions.txt
		for IND in $(cat ${TMPDIR}/01_download/failednames); do
		 	# Also remove data for these individuals we are re-running
		 	if [ -f /data/bioinfo/scratch/breakseq_fastqs/${DATE}_${NAME}/${IND}/selected_regions.fastq ] ; then
		 		rm /data/bioinfo/scratch/breakseq_fastqs/${DATE}_${NAME}/${IND}/selected_regions.fastq
		 	fi
			if [ -f ${TMPDIR}/01_download/readscount/${IND}.txt ] ; then
			 	rm ${TMPDIR}/01_download/readscount/${IND}.txt
			fi
			if [ -f ${LOGDIR}/01_download/${IND} ] ; then
			 	rm ${LOGDIR}/01_download/${IND}
			fi
		 	# Fill first failedregions
			for INV in $(cat $INVSFILE); do
				echo "${IND},${INV}" >> ${TMPDIR}/01_download/failedregions.txt
		 	done
		 	echo "${IND},Unmapped" >> ${TMPDIR}/01_download/failedregions.txt # Unmapped is also included, if it fails for n loops, it will show as failedregion, if it doesn't have reads it will take a little while, if it has reads it will be downloaded only once
		done

		# make recipe
		echo "executable = /bin/singularity
args = \"exec \\
--bind /data/bioinfo/scratch/breakseq_fastqs:/nfs/pic.es/user/r/rgomez/20210325_breakseq/${OUTDIR}/01_download/:rw \\
--bind /data/bioinfo/scratch/breakseq_bam:/nfs/pic.es/user/r/rgomez/20210325_breakseq/${DATADIR}/bamFiles/:rw \\
/data/bioinfo/software/rgomez_breakseq.sif \\
bash 20210325_breakseq/runall.sh download -r ${TMPDIR}/01_download/failedcases -s '\$(item)' -t '\$(item)' -n ${DATE}_${NAME}\"

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

			# MAKE LIST OF CASES FROM FAILEDREGIONS, and empty failedregions
			# The failedcases will be used to download only the necessary; the failedregions is filled in condor_submit with new information
			cat ${TMPDIR}/01_download/failednames > ${TMPDIR}/01_download/failednames_complete
			if [[ $MAXJOBS != "" ]]; then
				# head -n$MAXJOBS ${TMPDIR}/01_download/failednames_complete >  ${TMPDIR}/01_download/failednames
				# Make failednames with as many rows as MAXJOBS, since it is prepared to make a loop over samples, it whould work!
				NUMNAMES=$(wc -l ${TMPDIR}/01_download/failednames_complete | cut -d " " -f1)
				NUMCOLS=$(( ($NUMNAMES + $MAXJOBS -1) / $MAXJOBS)) # This is to make sure that it is rounded upwards
				awk '{ ORS = " " } { print }' ${TMPDIR}/01_download/failednames_complete | awk -v cols=$NUMCOLS '{for(i = cols+1; i<=NF; i+=cols) $i=RS $i}1' > ${TMPDIR}/01_download/failednames
			else
				cat ${TMPDIR}/01_download/failednames_complete >  ${TMPDIR}/01_download/failednames
			fi
			grep -v -f ${TMPDIR}/01_download/failednames_complete ${TMPDIR}/01_download/failedregions.txt > ${TMPDIR}/01_download/failedcases_complete # previous failed regions (that will not be run now)
			grep -f ${TMPDIR}/01_download/failednames_complete ${TMPDIR}/01_download/failedregions.txt > ${TMPDIR}/01_download/failedcases # previous failed regions (to test now )
			>${TMPDIR}/01_download/failedregions.txt # new failed regions (now empty)
			
			condor_submit ${LOGDIR}/01_download.sub
		
			waitForCondor 5 

			# at the end failedregions.txt is restored
			cat ${TMPDIR}/01_download/failedcases_complete >> ${TMPDIR}/01_download/failedregions.txt
			# CHECK CORRECT DOWNLOAD
			checkDownload 
		done

		rm ${TMPDIR}/01_download/failednames ${TMPDIR}/01_download/successnames

		# Make summary counts and copy them in folders
		cat ${TMPDIR}/01_download/readscount/* > ${OUTDIR}/01_download/readscount.txt
		cp ${OUTDIR}/01_download/readscount.txt /data/bioinfo/scratch/breakseq_fastqs/${DATE}_${NAME}_readscount.txt
		
	fi
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
	
	# import index path
	source ${DATADIR}/datos_librerias/specs.sh

	# Make analysis list
	if [[ $RESUME_B == "y" ]]; then
		echo "Resuming breakeq!"
		# CHECK CORRECT BREAKSEQ
		# -----------------------------------------------------

		# how many inds were downloaded?
		ls -l ${OUTDIR}/02_breakseq/*/*ref.sam | awk '$5>0{print $9}' | sed "s/^.*02_breakseq\///g" | sed 's/\/.*//g' > ${TMPDIR}/02_breakseq/successnames

		# Make a list of fails
		# This is the list for the queue
		grep -v -f ${TMPDIR}/02_breakseq/successnames $SAMPLESLIST  > ${TMPDIR}/02_breakseq/failednames
		# Counter
		TOTEST=$(cat ${TMPDIR}/02_breakseq/failednames | wc -l)


	else
		echo "Running breakseq!"
		cp $SAMPLESLIST ${TMPDIR}/02_breakseq/failednames
		TOTEST=1
	fi

	# set loop
	LOOPCOUNT=0
	MEMORY=2048
	
	while [[ $TOTEST -gt 0 ]]; do
		LOOPCOUNT=$(($LOOPCOUNT+1))
		echo "Starting loop $LOOPCOUNT"

	echo "executable = /bin/singularity
args = \"exec \\
--bind ${READ_FASTQS}:/nfs/pic.es/user/r/rgomez/20210325_breakseq/${OUTDIR}/01_download/:rw \\
--bind /data/bioinfo/common/${ASSEMBLY}:/nfs/pic.es/user/r/rgomez/20210325_breakseq/data/use/bowtie_index \\
--bind /data/bioinfo/scratch/breakseq_tmp/:/nfs/pic.es/user/r/rgomez/20210325_breakseq/${TMPDIR}/:rw \\
/data/bioinfo/software/rgomez_breakseq.sif \\
bash 20210325_breakseq/runall.sh breakseq -k $SCORE_MIN -m $COV_AROUND -l $LIBRARY_LENGTH -f ${OUTDIR}/01_download/ -s \$(Item) -t \$(Item) -n ${DATE}_${NAME}\"

output = ${LOGDIR}/02_breakseq/condor.out
error = ${LOGDIR}/02_breakseq/condor.err
log = ${LOGDIR}/02_breakseq/condor.log

request_cpus = 1
request_memory = $MEMORY


queue 1 from ${TMPDIR}/02_breakseq/failednames "  >  ${LOGDIR}/02_breakseq.sub
	

		condor_submit ${LOGDIR}/02_breakseq.sub
	
		waitForCondor 15 

		# CHECK CORRECT BREAKSEQ
		# -----------------------------------------------------

		# how many inds were downloaded?
		ls -l ${OUTDIR}/02_breakseq/*/*ref.sam | awk '$5>0{print $9}' | sed "s/^.*02_breakseq\///g" | sed 's/\/.*//g' > ${TMPDIR}/02_breakseq/successnames

		# Make a list of fails
		# This is the list for the queue
		grep -v -f ${TMPDIR}/02_breakseq/successnames $SAMPLESLIST  > ${TMPDIR}/02_breakseq/failednames

		# Counter
		TOTEST=$(cat ${TMPDIR}/02_breakseq/failednames | wc -l)

		#Memory
		echo "Duplicating memory to $(($MEMORY*2))"
		MEMORY=$(($MEMORY*2))



	done

	rm ${TMPDIR}/02_breakseq/failednames ${TMPDIR}/02_breakseq/successnames
	cat  ${OUTDIR}/02_breakseq/*/inisam_summary > ${OUTDIR}/02_breakseq/inisam_summary
	cat  ${OUTDIR}/02_breakseq/*/filsam_summary > ${OUTDIR}/02_breakseq/filsam_summary

	rm ${OUTDIR}/02_breakseq/*/inisam_summary ${OUTDIR}/02_breakseq/*/filsam_summary


	if [[ $KEEP_DOWNLOADS == "n" ]]; then
		# Delete used files from path
		echo "Deleting fastq files"
		ls /data/bioinfo/scratch/breakseq_fastqs/${DATE}_${NAME}/*/selected_regions.fastq | grep -f $SAMPLESLIST | rm
	fi
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
	# import index path
	source ${DATADIR}/datos_librerias/specs.sh

	echo "#!/bin/bash
INREF='${INREF}'
NOTINREF='${NOTINREF}'" > $DATADIR/03_processaligned_parameters.sh
	
	echo "executable = /bin/singularity
args = \"exec \\
--bind /data/bioinfo/common/${ASSEMBLY}:/nfs/pic.es/user/r/rgomez/20210325_breakseq/data/use/bowtie_index \\
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
		
	waitForCondor 1 

fi



#######################################################


echo "# Step 05 - TagSNPs
# =========================================================================== #"

if [[ $OVERRIDE -gt 5 ]] ; then
	echo "This step will not run now"
else

	mkdir -p ${LOGDIR}/05_tagsnps/ ${OUTDIR}/05_tagsnps/ ${TMPDIR}/05_tagsnps/
	
	# rm -r /data/bioinfo/scratch/breakseq_tmp/05_tagsnps/*
	# import index path
	source ${DATADIR}/datos_librerias/specs.sh

	echo "executable = /bin/singularity
args = \"exec \\
--bind /data/bioinfo/common/${ASSEMBLY}:/nfs/pic.es/user/r/rgomez/20210325_breakseq/data/use/bowtie_index:rw \\
--bind /data/bioinfo/scratch/breakseq_tmp/:/nfs/pic.es/user/r/rgomez/20210325_breakseq/${TMPDIR}/:rw \\
/data/bioinfo/software/rgomez_plots.sif \\
bash 20210325_breakseq/runall.sh tagsnps -m $MAX_ERROR -g $GENOTYPES_DIR -r \$(Item) -n ${DATE}_${NAME} -t \$(Item)  \"

output = ${LOGDIR}/05_tagsnps/condor.out
error = ${LOGDIR}/05_tagsnps/condor.err
log = ${LOGDIR}/05_tagsnps/condor.log

request_cpus = 1

queue 1 from ${INVSFILE}" >  ${LOGDIR}/05_tagsnps.sub

	condor_submit ${LOGDIR}/05_tagsnps.sub
		
	waitForCondor 1 

	echo -e "INV\tFILTER_genotypes\tFILTER_LD\tCHR\tPOS\tID\tLD"   > ${OUTDIR}/05_tagsnps/tagSNPs_gt08.txt
	cat /data/bioinfo/scratch/breakseq_tmp/05_tagsnps/*/tagSNPs_gt08.txt >> ${OUTDIR}/05_tagsnps/tagSNPs_gt08.txt
    echo -e "INV\tTYPE\tN_SAMPLES_${MAX_ERROR}\tchr_max_SNP_${MAX_ERROR}\tpos_max_SNP\tname_max_SNP\tLD_max_SNP" > ${OUTDIR}/05_tagsnps/tagSNPs_max.txt
  	cat /data/bioinfo/scratch/breakseq_tmp/05_tagsnps/*/tagSNPs_max.txt >> $OUTDIR/05_tagsnps/tagSNPs_max.txt

mkdir ${LOGDIR}/05_tagsnps_plot

echo "executable = /bin/singularity
args = \"exec \\
--bind /data/bioinfo/common/${ASSEMBLY}:/nfs/pic.es/user/r/rgomez/20210325_breakseq/data/use/bowtie_index:rw \\
--bind /data/bioinfo/scratch/breakseq_tmp/:/nfs/pic.es/user/r/rgomez/20210325_breakseq/${TMPDIR}/:rw \\
/data/bioinfo/software/rgomez_plots.sif \\
Rscript 20210325_breakseq/code/rscript/05_tagsnps.R 20210325_breakseq/${OUTDIR}/05_tagsnps/tagSNPs_max.txt 20210325_breakseq/${OUTDIR}/05_tagsnps/  \"

output = ${LOGDIR}/05_tagsnps_plot/condor.out
error = ${LOGDIR}/05_tagsnps_plot/condor.err
log = ${LOGDIR}/05_tagsnps_plot/condor.log

request_cpus = 1

queue 1 " >  ${LOGDIR}/05_tagsnps_plot.sub

	condor_submit ${LOGDIR}/05_tagsnps_plot.sub
		
	waitForCondor 1 
	

fi



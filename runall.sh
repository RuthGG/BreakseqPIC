#!/bin/bash
# Ruth Gómez Graciani
# 01 04 2020

###############################################################################
# Description:                                                                
# Run all the analyses in this project                                        
###############################################################################

# SAVE HISTORY pt 1
# Save date and command 
# =========================================================================== #
HISTORY="$@"

# USAGE 
# Help message. Don't forget to mention new options in README.md!!  
# =========================================================================== #

usage()
{
  echo "Usage: 

  $(basename $0) <COMMAND> [OPTIONS]
  "
  echo "Commands:
  download  -s -r -t -b -d -m -n -l      Download a set of reads from 1KGP samples, optionally run breakseq.
  breakseq  -f -s -t -m -n -l            Run breakseq on a list of fastqs.
  processaligned -s -a -t -m -n          Process reads that were aligned to probes.
  qualityanalysis -g -t -n               Make quality analysis for genotypes, including tagSNPs.
  "
  echo "Options:

    -h                                Show help.
    -s                                Samples list file.
    -r                                Regions list file (IDs for regions, e.g. inversion names). 
    -t                                Thread number.
    -b [y|n]                          Run breakseq now? (default n)
    -d [y|n]                          Delete tmp files now? (default n)
    -f                                Directory with fastq files (one directory per sample).
    -a                                Directory with aligned reads (breakseq results).
    -g                                Genotypes file or directory
    -m                                Parameter MIN. min. coverage for breakseq, min.alignment length in processaligned. Default 20_30, min. reads in quality control
    -n                                Optional parameter with process name. Default 'main'
    -l                                Library length (LIB_LEN parameter)
    -k                                Min score in bowtie2 for breakseq

 "

}


# PARSE OPTIONS
# Parse options and get help 
# =========================================================================== #

# Parse help message
while getopts "h" OPT; do
   case "$OPT" in
      h)  usage; exit 1 ;;
   esac
done
shift $((OPTIND -1))

# Save command, if any
COMMAND=$1; shift

# Set default optional variables
# CURDIR=$(pwd)     # Current directory
SCRIPTPATH="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
cd $SCRIPTPATH
STEP=00           # Steps

# Set empty mandatory variables
SAMPLES_FILE=""
REGIONS_FILE=""
THREAD="all"
B_OPTION="n"
D_OPTION="n"
FASTQ_DIR=""
ALIGNED_DIR=""
GENOTYPES_FILE=""
GENOTYPES_DIR=""
NAME="main"

MIN=30
MIN_COVER=20
LIB_LEN=300
# KGPTAG="n"
FILTER=0.03
# bowtie2 default
SCORE_MIN=L,-0.6,-0.6

# Parse command optons
case "$COMMAND" in
  #Download reads, optionally run breakeq
  download ) 
    while getopts "s:r:t:b:d:m:n:l:k:" OPTIONS ; do
      case "$OPTIONS" in
        s)  SAMPLES_FILE=${OPTARG} ;;
        r)  REGIONS_FILE=${OPTARG} ;;
        t)  THREAD=${OPTARG};;
        b)  B_OPTION=${OPTARG};;
        m)  MIN_COVER=${OPTARG};;
        l)  LIB_LEN=${OPTARG};;
        d)  D_OPTION=${OPTARG};;
        n)  NAME=${OPTARG};;
        k)  SCORE_MIN=${OPTARG};;
      esac
    done
    shift $((OPTIND -1))
  ;;
  # Run breakseq
  breakseq ) 
    while getopts "f:s:t:m:n:l:k:" OPTIONS ; do
      case "$OPTIONS" in
        s)  SAMPLES_FILE=${OPTARG} ;;
        t)  THREAD=${OPTARG};;
        m)  MIN_COVER=${OPTARG};;
        f)  FASTQ_DIR=${OPTARG};;
        l)  LIB_LEN=${OPTARG};;
        n)  NAME=${OPTARG};;
        k)  SCORE_MIN=${OPTARG};;
      esac
    done
    shift $((OPTIND -1))
  ;;
  # Process aligned reads (breakseq results)
  processaligned ) 
    while getopts "a:s:t:m:n:" OPTIONS ; do
      case "$OPTIONS" in
        s)  SAMPLES_FILE=${OPTARG} ;;
        t)  THREAD=${OPTARG};;
        a)  ALIGNED_DIR=${OPTARG};;
        m)  MIN=${OPTARG};;
        n)  NAME=${OPTARG};;
      esac
    done
    shift $((OPTIND -1))
  ;;
    # Make quality analysis for genotypes.
  qualityanalysis ) 
    while getopts "g:t:n:r:s:m:" OPTIONS ; do
      case "$OPTIONS" in
        g)  GENOTYPES_DIR=${OPTARG} ;;
        t)  THREAD=${OPTARG};;
        n)  NAME=${OPTARG};;
        r)  REGIONS_FILE=${OPTARG} ;;
        s)  SAMPLES_FILE=${OPTARG} ;;
        m)  FILTER=${OPTARG};;
      esac
    done
    shift $((OPTIND -1))
  ;;
esac


# Check that empty mandatory variables are full
if [ "$COMMAND" == "download" ]; then
  if [ -z "$SAMPLES_FILE" ] && [ -z "$REGIONS_FILE" ] ; then 
    echo "Remember that to use the '${COMMAND}' command, mandatory options are: -s, -r"; usage; exit
  fi
elif [ "$COMMAND" == "breakseq" ]; then
  if [ -z "$SAMPLES_FILE" ] && [ -z "$FASTQ_DIR" ] ; then 
    echo "Remember that to use the '${COMMAND}' command, mandatory options are: -s, -f"; usage; exit
  fi
elif [ "$COMMAND" == "processaligned" ]; then
  if [ -z "$SAMPLES_FILE" ] && [ -z "$ALIGNED_DIR" ] && [ -z "$MIN" ]; then 
    echo "Remember that to use the '${COMMAND}' command, mandatory options are: -s, -a, -m"; usage; exit
  fi
elif [ "$COMMAND" == "qualityanalysis" ]; then
  if [ -z "$GENOTYPES_DIR" ] && [ -z "$REGIONS_FILE" ] && [ -z "$SAMPLES_FILE" ]  ; then 
    echo "Remember that to use the '${COMMAND}' command, mandatory options are: -g, -r, -s"; usage; exit
  fi
else 
  usage; exit
fi

# SAVE HISTORY pt 2
# Save date and command 
# =========================================================================== #
DATE=$(date +%F)


# SAVE LOG pt 2
# Save date and command 
# =========================================================================== #
  exec 3>&1 4>&2
  trap 'exec 2>&4 1>&3' 0 1 2 3

# DOWNLOAD READS DATA
# =========================================================================== #
STEP=$(printf "%02d" $((${STEP}+1)))

if [ "$COMMAND" == "download" ]; then

  exec 1> analysis/${NAME}/log/${STEP}_${COMMAND}/${THREAD} 2>&1

  # Make directories
  TMPDIR="tmp/${NAME}/${STEP}_${COMMAND}"
  OUTDIR="analysis/${NAME}/${STEP}_${COMMAND}"
  DATADIR="analysis/${NAME}/data/"

  mkdir -p $OUTDIR $TMPDIR/readscount/

  # Take samples
  SAMPLES=$SAMPLES_FILE

  # Take regions
  echo "Take regions"
  # How to create a file with all invs
  REGIONS=$(cut -f1 $REGIONS_FILE)

  ## Loop per sample - align reads
  echo "Start looping samples"
  for SAMPLE in $SAMPLES; do
  
    echo "################"
    echo "$SAMPLE"
    echo "################"
    date
    # Make a folder for the sample
    mkdir -p ${OUTDIR}/${NAME}/${SAMPLE} 
        
    cd ${OUTDIR}/${NAME}/${SAMPLE} 

    # Remove output files if they exist
    if [ -f selected_regions.fastq ] ; then
      rm selected_regions.fastq
    fi
    if [ -f ${SCRIPTPATH}/${TMPDIR}/readscount/${SAMPLE}.txt ] ; then
      rm ${SCRIPTPATH}/${TMPDIR}/readscount/${SAMPLE}.txt
    fi
     

    # Loop per inversion - align reads
    for REGION in $REGIONS; do

      echo "---------------------"
      echo "$REGION Mapped reads "
      echo "---------------------"
      date
  

      # Locate bam file link
      # POP=$(grep $SAMPLE data/raw/1KGP_data/integrated_call_samples_v3.20130502.ALL.panel.txt | cut -f2)
      BAM_FILE=$(grep $SAMPLE ${SCRIPTPATH}/data/raw/1KGP_data/20130502.phase3.low_coverage.alignment.index | grep '\.mapped' | cut -f1 | sed -e 's/data.*alignment\///g')
      BAM_FILE="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/$SAMPLE/alignment/"$BAM_FILE

      # Download reads - Importante! Rango de extracción alrededor de la inversión! He puesto 20kb, pero quizá podría ser otro rango
      # Idealmente, habria que hacer una lista de coordenadas del alelo invertido para afinar esta region, pero todavia no existe
      CHR_REGION=$(grep "$REGION" ${SCRIPTPATH}/${DATADIR}/datos_librerias/bplib.coords | cut -f2 | uniq)
      START_REGION=$(grep  "$REGION" ${SCRIPTPATH}/${DATADIR}/datos_librerias/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | head -n+1)
      END_REGION=$(grep  "$REGION" ${SCRIPTPATH}/${DATADIR}/datos_librerias/bplib.coords | cut -f3,4 | sed -e 's/\t/\n/g' | sort | tail -n-1)
      echo "Original coords: "$CHR_REGION":"$START_REGION"-"$END_REGION
      
      START_REGION=$(($START_REGION-20000))
      END_REGION=$(($END_REGION+20000))
      echo "Sampled coords: "$CHR_REGION":"$START_REGION"-"$END_REGION
      
      # This includes a loop to resume download in case it was interrupted
      i=0
      l=0

      if [[ -z "$CHR_REGION" ]]; then
        echo "Empty coordinates"
        echo  "$SAMPLE","$REGION" >> ${SCRIPTPATH}/${TMPDIR}/failedregions.txt
      else

        while [ $i -eq 0 ]; do
          echo "# LOOP $l"

          if [[ $l -eq 1000 ]]; then
            echo "Maximum loop count! Aborting..."
            echo  "$SAMPLE","$REGION" >> ${SCRIPTPATH}/${TMPDIR}/failedregions.txt
          else

            l=$((l+1))
              

            ERRS=$(( samtools view $BAM_FILE $CHR_REGION":"$START_REGION"-"$END_REGION > tmp_download.txt ) 2>&1 )
            echo $ERRS

            # If tmp_download is not empty OR if tmp_download is empty but there were no errors and we tried more than 10 times already
            if [ -s tmp_download.txt ] || ( [ -z "$ERRS" ] && [ ! -s tmp_download.txt ] && [ $l -gt 10 ] ); then 
              echo "# Sucess!"
              cat tmp_download.txt | awk -v FS="\t" '{print "@" $1 "\n" $10 "\n+\n" $11}' >> selected_regions.fastq
              READS=$(($(cat tmp_download.txt | awk -v FS="\t" '{print "@" $1 "\n" $10 "\n+\n" $11}' | wc -l ) / 4)) 
              echo "Reads: " $READS
              echo "$SAMPLE,$REGION,$READS" >> ${SCRIPTPATH}/${TMPDIR}/readscount/${SAMPLE}.txt 
              i=1
            else
              echo "# Fail!"
              sleep 1
            fi

            echo "# Delete tmp file"
            rm tmp_download.txt
          fi
        done
      fi
    done
      # Unmapped reads

      echo "---------------------"
      echo " Unmapped reads "
      echo "---------------------"
      date
      # Locate bam file link
      # POP=$(grep $SAMPLE data/raw/1KGP_data/integrated_call_samples_v3.20130502.ALL.panel.txt | cut -f2)
      BAM_FILE=$(grep $SAMPLE ${SCRIPTPATH}/data/raw/1KGP_data/20130502.phase3.low_coverage.alignment.index | grep '\.unmapped' | cut -f1 | sed -e 's/data.*alignment\///g')
      BAM_FILE="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/phase3/data/$SAMPLE/alignment/"$BAM_FILE

      # This includes a loop to resume download in case it was interrupted
      i=0
      l=0
      while [ $i -eq 0 ]; do
        echo "# LOOP $l"
        l=$((l+1))

        # Work into tmp_download
        samtools view $BAM_FILE > tmp_download.txt
    
        if [ -s tmp_download.txt ]; then 
          # If file not empty, interrupt loop and concat to output file
          echo "# Success!"
          cat tmp_download.txt | awk -v FS="\t" '{print "@" $1 "\n" $10 "\n+\n" $11}' >> selected_regions.fastq
          READS=$(($(cat tmp_download.txt | awk -v FS="\t" '{print "@" $1 "\n" $10 "\n+\n" $11}' | wc -l ) / 4)) 
          echo "Reads: " $READS
          echo "$SAMPLE,Unmapped,$READS" >> ${SCRIPTPATH}/${TMPDIR}/readscount/${SAMPLE}.txt 
          i=1
        else
          echo "# Fail!"
          sleep 1
        fi

        echo "# Delete tmp file"
        rm tmp_download.txt

      done

      # Remove everything but the selected regions and return
      ls | grep -v 'selected_regions.fastq' | xargs rm 
      cd ${SCRIPTPATH}

      #  OBSOLETE! CHECK DIRECTORIES!!
      if [ "$B_OPTION" == "y" ]; then
        ## RUN BREAKSEQ
      
        echo "# Running BreakSeq on $SAMPLE."
            
        STEP_B=$(printf "%02d" $((${STEP}+1)))
        COMMAND_B="breakseq"

        OUTDIR_B="analysis/${STEP_B}_${COMMAND_B}/${DATE}/${NAME}"
        TMPDIR_B="tmp/${STEP_B}_${COMMAND_B}/${DATE}/${NAME}"

        echo "MIN_COVER = ${MIN_COVER}; LIB_LEN = ${LIB_LEN}" > ${OUTDIR_B}/README.txt

        mkdir -p project/logfiles/${STEP_B}_${COMMAND_B}_${DATE}_${NAME}/
        
        # Breakseq process 
        # TMPDIR=$1     
        # OUTDIR=$2  
        # FASTQ_DIR=$3   
        # MIN_COVER=$4
        # LIB_LEN=$5
        # SCORE_MIN=$6
        # SAMPLE=$7
        sh code/bash/02_runBreakseq.sh $TMPDIR_B $OUTDIR_B $FASTQ_DIR $MIN_COVER $LIB_LEN $SCORE_MIN $SAMPLE >> project/logfiles/${STEP_B}_${COMMAND_B}_${DATE}_${NAME}/${NAME}_${THREAD}


      fi
      
      # Remove and Return
      if [ "$D_OPTION" == "y" ]; then
        rm -r ${TMPDIR}/${SAMPLE}
      fi      
  done

fi


# BREAKSEQ
# =========================================================================== #
STEP=$(printf "%02d" $((${STEP}+1)))

if [ "$COMMAND" == "breakseq" ]; then

  exec 1> analysis/${NAME}/log/${STEP}_${COMMAND}/${THREAD} 2>&1

  # Make directories
  TMPDIR="tmp/${NAME}/${STEP}_${COMMAND}"
  OUTDIR="analysis/${NAME}/${STEP}_${COMMAND}"

  mkdir -p $OUTDIR $TMPDIR

  # Take samples
  SAMPLES=$SAMPLES_FILE

  ## Loop per sample - align reads
  echo "Start looping samples"
  for SAMPLE in $SAMPLES; do
      # TMPDIR=$1     
      # OUTDIR=$2  
      # FASTQ_DIR=$3   
      # MIN_COVER=$4
      # LIB_LEN=$5
      # SCORE_MIN=$6
      # SAMPLE=$7
      sh code/bash/02_runBreakseq.sh $TMPDIR $OUTDIR $FASTQ_DIR $MIN_COVER $LIB_LEN $SCORE_MIN $SAMPLE
  done

fi

# PROCESS ALIGNED READS
# =========================================================================== #
STEP=$(printf "%02d" $((${STEP}+1)))

if [ "$COMMAND" == "processaligned" ]; then

  exec 1> analysis/${NAME}/log/${STEP}_${COMMAND}/${THREAD} 2>&1

  # Make directories
  TMPDIR="tmp/${NAME}/${STEP}_${COMMAND}"
  OUTDIR="analysis/${NAME}/${STEP}_${COMMAND}"

  mkdir -p $OUTDIR $TMPDIR

  # Set variables
  REF_FASTA=data/use/bowtie_index/human_g1k_v37.fasta # Reference fasta
 # Library data before bowtie build
  DATADIR=analysis/${NAME}/data/datos_librerias

  source analysis/${NAME}/data/03_processaligned_parameters.sh

  # Take samples
  echo "Take samples"
  SAMPLES=$(cut -f1 $SAMPLES_FILE)
    
  # Indexamos la libreria en fasta
  if [ ! -e $DATADIR/bplib.fa.fai ]; then
          date
          echo "Indexing bplib.fa."
          echo
          samtools faidx $DATADIR/bplib.fa
  fi

  ## Loop per sample 
  echo "Start looping samples"
  for SAMPLE in $SAMPLES; do
  
    echo "################"
    echo "$SAMPLE"
    echo "################"

    # Make a folder for the sample
    mkdir -p ${TMPDIR}/${SAMPLE}
    
    echo "## STEP 1 MAKE BAM"
    # Archivos que contienen el número de reads que mapan de manera específica sobre
    # los puntos de rotura ('.uni') y de manera no específica ('.xun').
    # Para contar estos últimos como evidencia de la presencia del alelo estándard, debemos asegurarnos de que no son
    # reads repetitivos y que no mapan en ningún otro lugar. 
  
    # Pasamos el template del header a la carpeta del individuo concreto
      echo "  Creating header for $SAMPLE." 
      awk -v ID=$SAMPLE '{gsub(/XXX/, ID); print}' $DATADIR/ref.header.template > ${TMPDIR}/${SAMPLE}/ref.header

    # Tipos de archivos
    # .ini.sam -> archivo inicial con todos los alineamientos del breakseq sobre las librerias (INV y STD, o REF y DEL)
    # .fil.sam -> archivo filtrado del anterior con el minimo cover alreadedor del breakpoint
    # .ref.sam -> hits anteriores alineados contra el genoma de referencia (STD o REF)
    # .uni.sam -> hits finales unicos (que alinean en el invertido o el alelo alternativo al genoma de referencia)
    # .xun.sam -> hits finales no unicos (que alinean en el std o alelo de referencia del genoma)

    # Unir ref genome header con los reads que mapean en nuestras librerias en STD o REF (y tambien en el genoma) y los ordenamos
      echo "  Collecting and sorting alignments to the reference from individual $SAMPLE."
      if [ -f $ALIGNED_DIR/$SAMPLE/*.ref.sam ]; then
        awk '(FILENAME ~ /ref.header$/){print}((FILENAME ~ /sam$/) && (/^[^@]/)){print}' ${TMPDIR}/${SAMPLE}/ref.header $ALIGNED_DIR/$SAMPLE/*.ref.sam > ${TMPDIR}/${SAMPLE}/ref.sam
        samtools view -Sb  ${TMPDIR}/${SAMPLE}/ref.sam > ${TMPDIR}/${SAMPLE}/ref.bam
        samtools sort -n ${TMPDIR}/${SAMPLE}/ref.bam -o ${TMPDIR}/${SAMPLE}/ref.sorted.bam
        samtools view ${TMPDIR}/${SAMPLE}/ref.sorted.bam > ${TMPDIR}/${SAMPLE}/ref.sorted.sam
        # rm ${TMPDIR}/${SAMPLE}/ref.sam ${TMPDIR}/${SAMPLE}/ref.sorted.bam ${TMPDIR}/${SAMPLE}/ref.bam
                  
      else
        echo "  There are no hits to the human reference genome for $SAMPLE."
      fi

    # Pasamos el header de las librearias a la carpeta del individuo concreto
        echo "  Creating more headers for $SAMPLE."
        awk -v ID=$SAMPLE '{gsub(/XXX/, ID); print}' $DATADIR/bplib.header.template.txt > ${TMPDIR}/${SAMPLE}/header

    # Unir ref header de las libs con los reads que mapean en nuestras librerias en INV o alelo alternativo y los ordenamos
    # ADEMAS -> FILTRAR ALELOS STD, YA QUE SOLO LOS STD SON DE REFERENCIA! NO PUEDE HABER STD EN ESTE ARCHIVO! SOLO PUEDE INV!
      if [ -f $ALIGNED_DIR/$SAMPLE/*.uni.sam ]; then
        echo "  Collecting and sorting unique alignments to bplib from individual $SAMPLE."
        awk -v RG=$SAMPLE -v OFS="\t" -v MIN=$MIN '(FILENAME ~ /header$/){print}((FILENAME ~ /sam$/) && (/^[^@]/) && (length($10) >= MIN)){
          gsub(/:A$/,"",$3)
          print $0 "\tRG:Z:" RG}' ${TMPDIR}/${SAMPLE}/header $ALIGNED_DIR/$SAMPLE/*.uni.sam > ${TMPDIR}/${SAMPLE}/uni.sam

        ## FILTER ALELOS STD!
        # grep -v -P "\tSTDHsInv" ${TMPDIR}/${SAMPLE}/uni.sam | grep -v -P "\tREFHsInv" > ${TMPDIR}/${SAMPLE}/uni_aux.sam
        REFMATCH="^(${INREF})"
        awk -v pat=$REFMATCH '$3 !~ pat {print $0}' ${TMPDIR}/${SAMPLE}/uni.sam  > ${TMPDIR}/${SAMPLE}/uni_aux.sam

        rm ${TMPDIR}/${SAMPLE}/uni.sam
        mv ${TMPDIR}/${SAMPLE}/uni_aux.sam ${TMPDIR}/${SAMPLE}/uni.sam

        samtools view -Sb ${TMPDIR}/${SAMPLE}/uni.sam > ${TMPDIR}/${SAMPLE}/uni.bam
        samtools sort ${TMPDIR}/${SAMPLE}/uni.bam -o ${TMPDIR}/${SAMPLE}/uni.sorted.bam
        # rm ${TMPDIR}/${SAMPLE}/uni.sam ${TMPDIR}/${SAMPLE}/uni.bam
        samtools index ${TMPDIR}/${SAMPLE}/uni.sorted.bam
      else
        echo "  There are no unique hits to breakpoints for $SAMPLE."
      fi

    # Reads in xun.sam files map on breakpoints and also on the reference genome. For them
    # to be excluded as evidence of an alternative allele, they don't need to map with high
    # accuracy, nor uniquely. That's why lastly we are not using quality filters in breakseq.
    # However, for them to count as positive evidence of the reference allele, they need to
    # map to the coordinates of the reference breakpoint (not just anywhere in the genome)
    # and with some minimum quality.

    ## 1ª parte BLOQUE -> Unir ref header de las libs con los reads que mapean en nuestras librerias y en el genoma
    ## -> FILTRAR ALELOS INV, YA QUE SOLO LOS STD SON DE REFERENCIA! NO PUEDE HABER INV EN ESTE ARCHIVO! SOLO PUEDE STD!
      if [ -f $ALIGNED_DIR/$SAMPLE/*.xun.sam ]; then
        echo "  Collecting, filtering, and sorting non-unique alignments to bplib from individual $SAMPLE."
        awk -v RG=$SAMPLE -v OFS="\t" -v MIN=$MIN '(FILENAME ~ /header$/){print}((FILENAME ~ /sam$/) && (/^[^@]/) && (length($10) >= MIN)){
          gsub(/:A$/, "", $3)
          print $0 "\tRG:Z:" RG}' ${TMPDIR}/${SAMPLE}/header $ALIGNED_DIR/$SAMPLE/*.xun.sam > ${TMPDIR}/${SAMPLE}/xun.sam

        ## FILTER ALELOS INV!
        # grep -v -P "\INVHsInv" ${TMPDIR}/${SAMPLE}/xun.sam | grep -v -P "\tINSHsInv" | grep -v -P "\tDELHsInv" | grep -v -P "\tINSHsInv" | grep -v -P "\tALTHsInv" | grep -v -P "\tANCHsInv"> ${TMPDIR}/${SAMPLE}/xun_aux.sam
        NOREFMATCH="^(${NOTINREF})"
        awk -v pat=$NOREFMATCH '$3 !~ pat {print $0}'  ${TMPDIR}/${SAMPLE}/xun.sam > ${TMPDIR}/${SAMPLE}/xun_aux.sam


        rm ${TMPDIR}/${SAMPLE}/xun.sam
        mv ${TMPDIR}/${SAMPLE}/xun_aux.sam ${TMPDIR}/${SAMPLE}/xun.sam

    # Below, I am not filtering reads that happen in the ref bam file more than once, because
    # I want to respect the paired-ends, and because ambiguously mapped reads may have already
    # been filtered by the required minimum mapping quality of 15 (MinQual).

        MINQUAL=15

    #  filtrar reads que mapean en librerias si mapean en otras zonas del genoma tambien

        perl -we 'use strict; use warnings;
                  open(OUT, ">$ARGV[0]/xun.filtered.sam") || die "I cannot open xun.filtered.sam.\n";
                  open(BP, $ARGV[1]) || die "I cannot open bplib.coords.\n";
                  my (%BPchr, %BPstart, %BPend);
                  my $MinQual = $ARGV[2];

                  while (<BP>) {
                    chomp;
                    my @line = split /\t/, $_;
                    $BPchr{$line[0]} = $line[1];
                    $BPstart{$line[0]} = $line[2];
                    $BPend{$line[0]} = $line[3];
                  }
                  close BP;
                  open(REF, "$ARGV[0]/ref.sorted.sam") || die "I cannot open ref.sorted.sam.\n";
                  #my (%Freq, %Chr, %Pos, %Exclude);
                  my (%Chr, %Pos, %Exclude);
                  my $num = 0;
                  while (<REF>) {
                    $num++;
                    my @line = split /\t/, $_;
                    #$Freq{$line[0]}++;
                    #if ( $Freq{$line[0]} > 1 ) { $Exclude{$line[0]} = 1 }
                    my $good = 0;
                    if ($line[4] >= $MinQual) {
                      for my $k (keys %BPchr) {
                        ##############################################
                        ## Siguiente linea, no entiendo el +10 (¿?) ##
                        ##############################################
                        #if (($line[3] >= $BPstart{$k} + 10) && ($line[3] <= $BPend{$k}) && ($line[2] eq $BPchr{$k})) {
                        if (($line[3] >= $BPstart{$k}) && ($line[3] <= $BPend{$k}) && ($line[2] eq $BPchr{$k})) {
                          $good = 1;
                        }
                      }
                    }
                    #if ($good == 0) { $Exclude{$line[0]} = 1 }
                    if ($good == 0) { 
                      $Exclude{$line[0]} = 1 
                    }
                  }
    
                  my $excluded = scalar keys %Exclude;
                  print "On individual $ARGV[0], out of $num reads, $excluded have been excluded.\n";
                  close REF;
                  open(XUN, "$ARGV[0]/xun.sam") || die "I cannot open xun.sam.\n";


                  #########################################################
                  ## CUIDADO CON ESTE FILTRO - VERSION ANTIGUA DE NOMBRE ##
                  #########################################################
                  #while (<XUN>) {
                  #        if (/^@/) {print OUT $_}
                  #        else {
                  #                my @line = split /\t/, $_;
                  #                unless ((exists $Exclude{$line[0]}) || ($line[2] =~ /^INV/) || !($line[2] =~ /REF/)) {print OUT $_}
                  #        }
                  #}
                  while (<XUN>) {
                    if (/^@/) {
                      print OUT $_ ;
                    }
                    else {
                      my @line = split /\t/, $_;
                      unless (exists $Exclude{$line[0]}) {
                        print OUT $_ ;
                      }
                    }
                  }
                  close XUN;' ${TMPDIR}/${SAMPLE}/ $DATADIR/bplib.coords $MINQUAL

    # Ordenar reads
        if grep -q -v -P "^@" ${TMPDIR}/${SAMPLE}/xun.filtered.sam; then
          samtools view -Sb ${TMPDIR}/${SAMPLE}/xun.filtered.sam > ${TMPDIR}/${SAMPLE}/xun.filtered.bam
          samtools sort ${TMPDIR}/${SAMPLE}/xun.filtered.bam -o ${TMPDIR}/${SAMPLE}/xun.sorted.bam
          rm ${TMPDIR}/${SAMPLE}/xun.filtered.bam
          samtools index ${TMPDIR}/${SAMPLE}/xun.sorted.bam
        else
          echo "  Indeed, it seems that all reads were excluded..."
        fi
        rm ${TMPDIR}/${SAMPLE}/xun.sam ${TMPDIR}/${SAMPLE}/xun.filtered.sam


      else
        echo "  There are no non-unique hits to the breakpoints for $SAMPLE."
      fi

    echo "## STEP 2 FILTERUNIXUN"
    # After processing of the mapped reads, the available evidence
    # to predict the genotypes is filtered, cleaned, and it is now that the genotypes
    # should be predicted, using uni.sorted.bam and xun.sorted.bam, if available.
      echo "  Filtering uni.sorted.bam"
      if [ -e ${TMPDIR}/${SAMPLE}/uni.sorted.bam ] ; then
          samtools view ${TMPDIR}/${SAMPLE}/uni.sorted.bam > ${TMPDIR}/${SAMPLE}/uni.sorted.sam
          python2 /breakseq-1.3/bin/svMap/svhit.py ${TMPDIR}/${SAMPLE}/uni.sorted.sam > ${TMPDIR}/${SAMPLE}/filtered.uni
          # rm ${TMPDIR}/${SAMPLE}/uni.sorted.sam
      fi
      echo "  Filtering xun.sorted.bam"
      if [ -e ${TMPDIR}/${SAMPLE}/xun.sorted.bam ] && [ ! -e ${TMPDIR}/${SAMPLE}/filtered.xun ]; then
          samtools view ${TMPDIR}/${SAMPLE}/xun.sorted.bam > ${TMPDIR}/${SAMPLE}/xun.sorted.sam
          python2 /breakseq-1.3/bin/svMap/svhit.py  ${TMPDIR}/${SAMPLE}/xun.sorted.sam > ${TMPDIR}/${SAMPLE}/filtered.xun
          # rm ${TMPDIR}/${SAMPLE}/xun.sorted.sam
      fi

      echo "  Run interpreting perl script"
      perl -we 'use strict; use warnings;

                my $iOpt= $ARGV[0];
                my $sOpt= $ARGV[1];

                #############################
                #### Print header
                #############################

                #print  "#Inversion\tAllele\tBreakpoint\tIndividual\tReads\n";
                # I do not print header, to be able to sort the result later.

                #################################################
                #### Open files .uni and .xun for each individual
                #################################################

                if (-e $iOpt . ".uni") {
                  open(UNI, ($iOpt . ".uni")) or die "I cannot open $iOpt.uni\n";
                  while (<UNI>) {
                          my @line = split /\t/, $_;
                          if ($line[0] =~ /^([A-Z]{3})(\D*0*\d+)BP(0|1|2|3|4|5|6|7|8|9)(.+)$/) {
                            print  $2, "\t", $1, "\tBP", $3,$4, "\t", $sOpt, "\t", $line[1], "\n";
                          }
                  }
                  close(UNI);
                }

                if (-e $iOpt . ".xun") {
                  open(XUN, ($iOpt . ".xun")) or die "I cannot open $iOpt.xun.\n";
                  while (<XUN>) {
                          my @line = split /\t/, $_;
                          if ($line[0] =~ /^([A-Z]{3})(\D*0*\d+)BP(0|1|2|3|4|5|6|7|8|9)(.+)$/) {
                            print  $2, "\t", $1, "\tBP", $3,$4, "\t", $sOpt, "\t", $line[1], "\n";
                          }
                            
                  }
                  close(XUN);
                };' ${TMPDIR}/${SAMPLE}/filtered $SAMPLE > ${TMPDIR}/${SAMPLE}/z1
      
  done

  # This is technically step 3.1

  cat ${TMPDIR}/*/z1 > $OUTDIR/Results_reads.txt


  Rscript code/rscript/03_inferGenotypes.R analysis/${NAME}/


  rm -r ${TMPDIR}


fi



# QUALITY ANALYSIS
# =========================================================================== #
STEP=$(printf "%02d" $((${STEP}+1)))

if [ "$COMMAND" == "qualityanalysis" ]; then
  
  exec 1> analysis/${NAME}/log/${STEP}_${COMMAND}/${THREAD} 2>&1

 # Make directories
  TMPDIR="tmp/${NAME}/${STEP}_${COMMAND}"
  OUTDIR="analysis/${NAME}/${STEP}_${COMMAND}"
  FASTA="analysis/${NAME}/data/datos_librerias/bplib.fa"

  mkdir -p $OUTDIR $TMPDIR

  
  # Variables not imputed
  # Based on a static directory structure
  REF_GENOTYPES="data/raw/GlobalInvGenotypes_v3.2_132Invs_20210528_Genotypes.csv"
  DATADIR="analysis/${NAME}/data/datos_librerias/"

  READS_FILE="${GENOTYPES_DIR}/Results_reads.txt"
  GENOTYPES_FILE="${GENOTYPES_DIR}/GTypes_FinalDataSet.txt"

  ############################################
  ## MAKE PROBE QUALITY CONTROL             ##
  ############################################

  cat analysis/${NAME}/02_breakseq/*/*uni.sam | cut -f3 > ${TMPDIR}/unisam_summary
  cat analysis/${NAME}/02_breakseq/*/*xun.sam | cut -f3> ${TMPDIR}/xunsam_summary
  cp analysis/${NAME}/02_breakseq/inisam_summary  ${TMPDIR}/inisam_summary
  cp analysis/${NAME}/02_breakseq/filsam_summary  ${TMPDIR}/filsam_summary

  Rscript code/rscript/04_filterSummary.R ${TMPDIR}/ $READS_FILE $OUTDIR $FASTA

  rm  ${TMPDIR}/unisam_summary ${TMPDIR}/xunsam_summary ${TMPDIR}/filsam_summary ${TMPDIR}/inisam_summary

  ############################################
  ## COMPARE BS GTYPES WITH INVFEST PROJECT ##
  ############################################
  # echo "Rscript code/rscript/04_qualityanalysis.R $GENOTYPES_FILE $SAMPLES_FILE $REF_GENOTYPES $OUTDIR $REGIONS_FILE"
  Rscript code/rscript/04_qualityanalysis.R $GENOTYPES_FILE $SAMPLES_FILE $REF_GENOTYPES $OUTDIR $REGIONS_FILE $FILTER

  GENOTYPES_FILE="${OUTDIR}/genotypesClean.txt"


  ########################
  ## CALCULATE TAG SNPS ##
  ########################
 
  bash code/bash/05_tagsnps.sh $TMPDIR $OUTDIR $REGIONS_FILE $FILTER $GENOTYPES_FILE $DATADIR $SCRIPTPATH
  Rscript code/rscript/05_tagsnps.R $OUTDIR/tagSNPs_max.txt $OUTDIR/

  rm $TMPDIR
fi


#!/bin/bash

# Set here all the variables for runCommands.sh

# Name of the process 
# ALERT!!! YOU WON'T BE WARNED IF THIS ALREADY EXISTS!
# =========================================================================== #

NAME="agdp"

# Manually set date, comment if we want it automatic
DATE="2022-04-18"


# To override some parts

OVERRIDE=2 # This is the first process number we want to run, e.g. 2 will not run 0, and 1

# To override:
	# 0 library -> nothing needed, library is always in the same place /data/bioinfo/scratch/breakseq_fastqs
	# 1 download -> we can override it by just having a FASTQ variable, even with 0 value override
		# If we want to override the download, we can set a directory with previous downloads
		# It has to be the directory with individual directories
		# READ_FASTQS="/data/bioinfo/scratch/breakseq_fastqs/2022-02-25_proves_v2.3.3"
	# 2 breakseq
		#BREAKSEQ_RESULTS="analysis/" 
	# 3 alignment analysis - it can be automatically assigned!
		# GENOTYPES_DIR="analysis/2021-06-30_basicRun/03_processaligned/"
	# 4 quality control and tagnsps
		# It is automatically assigned!



# STEP 00 - Library
# =========================================================================== #

# library seed name
# this name must exist in data/raw/seed_librerias
LIBRARY_SEED="v2.3.3.300"

# This can be used in the future to automatically trim the original probes
LIBRARY_LENGTH=300

# Do we want to update all library or just to make the base files?
LIBRARY_UPDATE="no"

# STEP 01 - Download
# =========================================================================== #

# individuals to test"
SAMPLESLIST="data/use/ids_agdp.txt"
SAMPLESFILE="data/use/agdp_data/pathIndex_agdp.txt"

# inversions to test - it can be empty to analyze all
INVSFILE=""

# Set whether I want to delete individual files (y or n)
KEEP_DOWNLOADS="y"

# Resume download with the same name (y or n)
# To resume download we need a TMP dir with a path equivalent to this process that contains a readscount summary for each individual
RESUME="n"


# Step 02 - Breakseq
# =========================================================================== #

# Coverage around the breakpoint that we want to take into account
# Later I can modify it to be a list
COV_AROUND=10

# --score-min option in BOWTIE2 for the breakseq run
SCORE_MIN=L,0,-0.3

# Step 03 - ProcessAligned
# =========================================================================== #

# Later I can modify it to be a list
MIN_LENGTH=20

# Allele names, separated with |
INREF="STD|REF"
NOTINREF="INV|INS|DEL|ALT|ANC"

# Step 04 - QualityAnalysis
# =========================================================================== #
# Later I can modify this to accept several breakseq + alignment results

# the maximum error accepted (<=)
MAX_ERROR=0.03


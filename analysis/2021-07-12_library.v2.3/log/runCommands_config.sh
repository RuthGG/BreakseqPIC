#!/bin/bash

# Set here all the variables for runCommands.sh

# Name of the process 
# ALERT!!! YOU WON'T BE WARNED IF THIS ALREADY EXISTS!
# =========================================================================== #

NAME="library.v2.3"

# Manually set date, comment if we want it automatic
DATE="2021-07-12"

# To override some parts

OVERRIDE=4 # This is the first process number we want to run, e.g. 2 will not run 0, and 1

# To override:
	# 0 library -> nothing needed, library is always in the same place /data/bioinfo/scratch/breakseq_fastqs
	# 1 download -> we can override it by just having a FASTQ variable, even with 0 value override
		# If we want to override the download, we can set a directory with previous downloads
		# It has to be the directory with individual directories
		READ_FASTQS="/data/bioinfo/scratch/breakseq_fastqs/01_download/2021-05-07/"
	# 2 breakseq
		# BREAKSEQ_RESULTS="analysis/02_breakseq/2021-06-03/smallinvs_min_20" 
	# 3 alignment analysis - it can be automatically assigned!
		# GENOTYPES_DIR="analysis/2021-06-30_basicRun/03_processaligned/"
	# 4 quality control and tagnsps
		# It is automatically assigned!



# STEP 00 - Library
# =========================================================================== #

# library seed name
# this name must exist in data/raw/seed_librerias
LIBRARY_SEED="v2.3.1.300"

# This can be used in the future to automatically trim the original probes
LIBRARY_LENGTH=300

# Do we want to update all library or just to make the base files?
LIBRARY_UPDATE="yes"

# STEP 01 - Download
# =========================================================================== #

# individuals to test
SAMPLESFILE="data/use/static_dataset_list.txt"

# inversions to test - it can be empty to analyze all
INVSFILE=""

# Set whether I want to delete individual files (y or n)
KEEP_DOWNLOADS="y"


# Step 02 - Breakseq
# =========================================================================== #

# Coverage around the breakpoint that we want to take into account
# Later I can modify it to be a list
COV_AROUND=20 

# Step 03 - ProcessAligned
# =========================================================================== #

# Later I can modify it to be a list
MIN_LENGTH=30

# Step 04 - QualityAnalysis
# =========================================================================== #
# Later I can modify this to accept several breakseq + alignment results

# The minimum amount of reads to filter each individual to be used
MIN_READS=10


# AÑADIR LISTA DE ACEPTADOS PARA UNI Y ACEPTADOS PARA XUN
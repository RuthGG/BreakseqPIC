###############################################################################
# Description:                                                                
# Commands to use in this project        
###############################################################################

# =========================================================================== #
# TO USE SINGULARITY ENVIRONMENTS 
# =========================================================================== #

#----------------------------------
# 1st, open interactive conda shell
#----------------------------------
condor_submit -interactive

#----------------------------------
# 2nd, create singularity comand
# ------------------------------

# Enter singularity shell OR singularity exec command 

	singularity shell \
	singularity exec \

# Add bind folders as needed --> origin:existing dir:[rw]

	# Bowtie indexes
	--bind /data/bioinfo/common/bowtie2_index:/nfs/pic.es/user/r/rgomez/20210325_breakseq/data/use/bowtie_index 

	# SMALL Bowtie indexes
	--bind /data/bioinfo/common/bowtie2_index_150:/nfs/pic.es/user/r/rgomez/20210325_breakseq/data/use/bowtie_index 

	# FASTQ files to download
	--bind /data/bioinfo/scratch/breakseq_fastqs:/nfs/pic.es/user/r/rgomez/20210325_breakseq/tmp:rw

	# FASTQ files to use
	--bind /data/bioinfo/scratch/breakseq_fastqs/01_download/2021-05-07:/nfs/pic.es/user/r/rgomez/20210325_breakseq/tmp/fastqs/:rw

	# 1KGP VCFs
	--bind /data/bioinfo/scratch/1KGP_hg19_raw/:/nfs/pic.es/user/r/rgomez/20210325_breakseq/data/raw/1KGP_data/vcf

# Add container name

	/data/bioinfo/software/rgomez_breakseq.sif 

# An example of a complete basic command
	singularity shell data/bioinfo/software/rgomez_breakseq.sif 

#------------------------------------------------------
# Alternatively, execute from conda, as explained below
# -----------------------------------------------------

# =========================================================================== #
# TO UPDATE SINGULARITY ENVIRONMENTS 
# =========================================================================== #

# Change in local
# Push

# As in:
# -----------------------------------------------------
	
	git add . ;git commit -m "<MESSAGE>"; git push

# After ~ 10 mins ( 428 seconds ):
# cd /data/bioinfo/software
# rm rgomez_breakseq.sif
# singularity pull --disable-cache shub://jmurga/bgd-pic:breakseq
singularity pull --disable-cache --arch amd64 library://ruthgg/default/breakseq:latest 

# As in:
# -----------------------------------------------------

	cd /data/bioinfo/software; rm rgomez_breakseq.sif; singularity pull --disable-cache shub://jmurga/bgd-pic:breakseq
	cd /data/bioinfo/software; rm rgomez_breakseq.sif;  singularity pull --disable-cache library://ruthgg/default/breakseq:sha256.1a436eaca59c7a5026e9461f6332ef4b0860abbfea84617ff347aac4046f4753; mv breakseq_sha256.1a436eaca59c7a5026e9461f6332ef4b0860abbfea84617ff347aac4046f4753.sif rgomez_breakeq.sif
	cd /data/bioinfo/software; rm rgomez_plots.sif;  
	singularity pull --disable-cache --arch amd64 library://ruthgg/default/plots:latest;mv plots_latest.sif rgomez_plots.sif 

# =========================================================================== #
# TO USE SRA TOOLS
# =========================================================================== #
singularity pull --disable-cache docker://ncbi/sra-tools
condor_submit -interactive
singularity shell \
--bind /data/bioinfo/scratch/breakseq_fastqs/SRAfiles/:/nfs/pic.es/user/r/rgomez/local-file-caching/:rw \
/data/bioinfo/software/sra-tools_latest.sif

# =========================================================================== #
# SETUP THE MAIN DIR
# =========================================================================== #

#------------------------------------------------------
# To make a list of individuals
# -----------------------------------------------------

# Split into a number of files (e.g. = 400)
	
	split --number=l/400 data/use/static_dataset_list.txt  data/use/sampleLists/samples_1kgp_ --numeric-suffixes=0 --suffix-length=1 --additional-suffix='.txt'

# Split with n number of lines per file (e.g. = 1) --> CURRENT!

	split -l 1 data/use/static_dataset_list.txt  data/use/sampleLists/samples_1kgp_ --numeric-suffixes=0 --suffix-length=3 --additional-suffix='.txt'

#------------------------------------------------------
# To download the 1KGP VCFs
# This is not exactly correct! ChrX works different, make conditional!!
# -----------------------------------------------------

  	for CHR in {1..22} X ;do
      echo $CHR
       wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz
       wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/ALL.chr${CHR}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi
    done

# =========================================================================== #
# To make a RUNALL.SUB file
# =========================================================================== #

#------------------------------------------------------
# TEMPLATE
# -----------------------------------------------------

executable = /bin/singularity
args = "exec \
/data/bioinfo/software/rgomez_breakseq.sif \
bash 20210325_breakseq/runall.sh "

output = condor.out
error = condor.err
log = condor.log

request_cpus = 1

queue 

#------------------------------------------------------
# QUEUE EXAMPLES
# -----------------------------------------------------

queue from seq -w 000 413 |
queue 1 in (110,123,148,205,219,373,413)


# =========================================================================== #
# RUNALL.SUB LAST VERSIONS
# To run this directly on an interactive condor shell, enter a singularity
# shell, or exec from interactive condor, replacing $Item manually.
# =========================================================================== #

#------------------------------------------------------
# ONLY DOWNLOAD - keeping files (default)
# NO NAME TESTED
# -----------------------------------------------------

executable = /bin/singularity
args = "exec \
--bind /data/bioinfo/scratch/breakseq_fastqs:/nfs/pic.es/user/r/rgomez/20210325_breakseq/tmp:rw \
/data/bioinfo/software/rgomez_breakseq.sif \
bash 20210325_breakseq/runall.sh download -r data/use/inversions_completeList.txt -s data/use/sampleLists/samples_1kgp_$(Item).txt -t $(Item)"

output = condor.out
error = condor.err
log = condor.log

request_cpus = 1

queue from seq -w 000 413 |

	# TO CHECK CORRECT DOWNLOAD
	# -----------------------------------------------------

	cat  data/use/static_dataset_list.txt  | cut -f1 > test
	# how many inds were downloaded? - two ways (check they are the same)
	ls /data/bioinfo/scratch/breakseq_fastqs/01_download/2021-05-07/ | wc -l 
	cat project/logfiles/2021-05-07_download/* > project/logfiles/2021-05-07_download/00_download_all
	grep 'NA[0-9]*$' project/logfiles/2021-05-07_download/00_download_all | wc -l
	# Take the list with less names
	# Make a list of fails
	[command with the less names] > test_b
	grep -v -f test_b test  > failed
	grep -n -f failed data/use/static_dataset_list.txt | sed 's/:.*//' | awk '{print ($1-1)}' | tr '\n' ','
	# This is the list for the queue
	# Also check that all the dirs contain the same files

#------------------------------------------------------
# DOWNLOAD + BREAKSEQ - keeping files, MIN_COV = 20
# can delete files if we add '-d y'
# NO NAME TESTED
# -----------------------------------------------------

executable = /bin/singularity
args = "exec \
--bind /data/bioinfo/scratch/breakseq_fastqs:/nfs/pic.es/user/r/rgomez/20210325_breakseq/tmp:rw \
/data/bioinfo/software/rgomez_breakseq.sif \
bash 20210325_breakseq/runall.sh download -b y -m 20 -r data/use/inversions_completeList.txt -s data/use/sampleLists/samples_1kgp_$(Item).txt -t $(Item)"

output = condor.out
error = condor.err
log = condor.log

request_cpus = 1

queue from seq -w 000 413 |

#------------------------------------------------------
# BREAKSEQ - MIN_COV = 20
# To update, remember to change the date from data download
# -----------------------------------------------------

executable = /bin/singularity
args = "exec \
--bind /data/bioinfo/scratch/breakseq_fastqs/01_download/2021-05-07:/nfs/pic.es/user/r/rgomez/20210325_breakseq/tmp/fastqs/:rw \
--bind /data/bioinfo/common/bowtie2_index:/nfs/pic.es/user/r/rgomez/20210325_breakseq/data/use/bowtie_index \
/data/bioinfo/software/rgomez_breakseq.sif \
bash 20210325_breakseq/runall.sh breakseq -m 20 -l 300 -f tmp/fastqs -s data/use/sampleLists/samples_1kgp_$(Item).txt -t $(Item) -n smallinvs_min_20"

output = condor.out
error = condor.err
log = condor.log

request_cpus = 1

queue from seq -w 000 413 |

#------------------------------------------------------
# BREAKSEQ - SHORT  INVS
# -----------------------------------------------------

executable = /bin/singularity
args = "exec \
--bind /data/bioinfo/scratch/breakseq_fastqs/01_download/2021-05-07:/nfs/pic.es/user/r/rgomez/20210325_breakseq/tmp/fastqs/:rw \
--bind /data/bioinfo/common/bowtie2_index_150:/nfs/pic.es/user/r/rgomez/20210325_breakseq/data/use/bowtie_index \
 /data/bioinfo/software/rgomez_breakseq.sif \
bash 20210325_breakseq/runall.sh breakseq -m 10 -l 150 -f tmp/fastqs -s data/use/sampleLists/samples_1kgp_$(Item).txt -t $(Item) -n smallinvs_min_10"

output = condor.out
error = condor.err
log = condor.log

request_cpus = 1

queue from seq -w 000 413 |

	# TO CHECK CORRECT BREAKSEQ
	# -----------------------------------------------------
	# Check if important file is 0
	ls -la  analysis/02_breakseq/2021-05-10/*/*uni.sam  | awk '$5 == 0{print $9}'| sed 's/analysis\/02_breakseq\/2021-05-10\///g' | sed 's/\/.*//g' > failed

	grep -n -f failed data/use/static_dataset_list.txt | sed 's/:.*//' | awk '{print ($1-1)}' | tr '\n' ','
	# This is the list for the queue

	ls analysis/02_breakseq/2021-05-25/smallinvs_min_20 | grep NA* > success
	grep -n -v -f success data/use/static_dataset_list.txt | sed 's/:.*//' | awk '{print ($1-1)}' | tr '\n' ','

# -----------------------------------------------------
# ALIGNED READS ANALYSIS - MIN = 30,70
# To update, remember to change the date from the breakseq
# -----------------------------------------------------

executable = /bin/singularity
args = "exec \
--bind /data/bioinfo/common/bowtie2_index:/nfs/pic.es/user/r/rgomez/20210325_breakseq/data/use/bowtie_index \
/data/bioinfo/software/rgomez_breakseq.sif \
bash 20210325_breakseq/runall.sh processaligned -a analysis/02_breakseq/2021-05-11/min_20/ -s data/use/static_dataset_list.txt -m $(Item) -t min_20_$(Item) -n min_20_$(Item)"

output = condor.out
error = condor.err
log = condor.log

request_cpus = 1

queue 1 in (30, 70)

# -----------------------------------------------------
# QUALITY CONTROL
# To update, remember to change the date from the alignment
# -----------------------------------------------------

executable = /bin/singularity
args = "exec \
--bind /data/bioinfo/scratch/1KGP_hg19_raw/:/nfs/pic.es/user/r/rgomez/20210325_breakseq/data/raw/1KGP_data/vcf \
/data/bioinfo/software/rgomez_plots.sif \
bash 20210325_breakseq/runall.sh qualityanalysis \
-g analysis/03_processaligned/2021-06-04/smallinvs_min_$(min1)_$(min2)/GTypes_FinalDataSet.txt -s data/use/static_dataset_list.txt -r data/use/small_invs.txt -m 10 -n smallinvs_min_$(min1)_$(min2) -t smallinvs_min_$(min1)_$(min2) "



output = condor.out
error = condor.err
log = condor.log

request_cpus = 1

queue min1, min2 from (
  10,10
  10,90
  10,30
  20,10
  20,90
  20,30
  30,10
  30,90
  30,30
)






# -------- scratch

du -h -d 1 project/logfiles/

singularity exec \
--bind /data/bioinfo/scratch/breakseq_fastqs/01_download/2021-05-07:/nfs/pic.es/user/r/rgomez/20210325_breakseq/tmp/fastqs/:rw \
--bind /data/bioinfo/common/bowtie2_index:/nfs/pic.es/user/r/rgomez/20210325_breakseq/data/use/bowtie_index \
/data/bioinfo/software/rgomez_breakseq.sif \
bash 20210325_breakseq/runall.sh breakseq -m 20 -f tmp/fastqs -s data/use/sampleLists/samples_1kgp_000.txt -t 000 -n min_20


# -----------------------------------------------------
# TAG SNP ANALYSIS - DEPRECATED _ FUSED WITH QUALITY CONTROL
# Remember to update date from breakseq results
# -----------------------------------------------------

executable = /bin/singularity
args ="exec \
 --bind /data/bioinfo/common/bowtie2_index:/nfs/pic.es/user/r/rgomez/20210325_breakseq/data/use/bowtie_index \
 --bind /data/bioinfo/scratch/1KGP_hg19_raw/:/nfs/pic.es/user/r/rgomez/20210325_breakseq/data/raw/1KGP_data/vcf \
 /data/bioinfo/software/rgomez_breakseq.sif \
bash 20210325_breakseq/runall.sh tagsnps -g analysis/03_processaligned/2021-05-12/min_20_30/GTypes_FinalDataSet.txt -r data/use/inversions_completeList.txt -n min_20_30"


output = condor.out
error = condor.err
log = condor.log

request_cpus = 1

queue 1 
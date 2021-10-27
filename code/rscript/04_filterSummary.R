rm(list=ls()) 
#####################
# Ruth Gómez Graciani
# 06 07 2021

###############################################################################
# Description:                                                                
# Make quality control to know which probes are not working         
###############################################################################

# LOAD ARGUMENTS
# =========================================================================== #
args = commandArgs(trailingOnly=TRUE)

# Example
# args[1]<-"tmp/2021-09-29_v2.3.2_deftest/02_breakseq/" # Path with ini fil uni xun
# args[2]<-"analysis/2021-09-29_v2.3.2_deftest/03_processaligned/Results_reads.txt"  # Final read counts
# args[3]<-"analysis/2021-09-29_v2.3.2_deftest/04_qualityanalysis/" #Output path
# args[4]<-"analysis/2021-09-29_v2.3.2_deftest/data/datos_librerias/bplib.fa"
# 


# # Test if there is at least one argument: if not, return an error
if (length(args)<4) {
  stop("One input file, one sample names file, one reference genotypes file, one output directory, a coordinates file", call.=FALSE)
}



# LOAD PACKAGES
# =========================================================================== #

# library(ggplot2)
# library(reshape2)
# library(grid)


#####################

pathname<-args[1]

tablelist<-list()

for (i in c("ini", "fil", "xun", "uni")) {
  file <- paste0(pathname, i, "sam_summary")
  filesam<-read.table(file, stringsAsFactors = F)
  
  if (is.null(filesam$V2)){
    filesam$V2 <-1
  }
  
  filesam<-aggregate( V2  ~ V1, filesam, sum)
  
  colnames(filesam)<-c("probe", paste0("count.",i))

  tablelist[[i]]<-filesam
}

fasta<-read.table(args[4], stringsAsFactors = F)
fasta<-fasta[grep("^>", fasta$V1),]
fasta<-as.data.frame(sub(">", "", fasta ))
colnames(fasta)<-"probe"
tablelist[["fasta"]]<-fasta

summary_probes<-Reduce(function(x,y) merge(x = x, y = y, by = "probe", all=TRUE), 
       tablelist)

# Now with reads
# summary_probes$code<-substr(summary_probes$probe, 1, 15)

reads<-read.table(args[2])
reads$code<-paste0(reads$V2, reads$V1, reads$V3)

reads_counted<-aggregate( V5  ~ code, reads, sum)
colnames(reads_counted)<-c("probe", "count.qc")

summary_probes<-merge(summary_probes, reads_counted, by = "probe", all = TRUE)

# Look for errors!
# summary_probes$code[duplicated(summary_probes$code)]

summary_probes[is.na(summary_probes)]<-0
summary_probes$code<-NULL


# Mark statistically divergent probes

# without probes 
summary_probes$diagnosis1<- NA
summary_probes[which((summary_probes$count.xun == 0) & (substr(summary_probes$probe, 1, 3) %in% c("STD", "REF"))), "diagnosis1"]<-"emptyRef"
summary_probes[which((summary_probes$count.uni == 0) & !(substr(summary_probes$probe, 1, 3) %in% c("STD", "REF"))), "diagnosis1"]<-"emptyAlt"

#  undetected or filtered
summary_probes$diagnosis2<-ifelse(summary_probes$count.ini == 0, "failedDetection", 
                                  ifelse(summary_probes$count.fil == 0, "allFiltered", 
                                         ifelse(summary_probes$count.qc == 0, "failedQC", NA)))

# Error match
summary_probes$diagnosis3<-ifelse( substr(summary_probes$probe, 1, 3) %in% c("STD", "REF"),
                                   # if is std,
                                   ifelse(summary_probes$count.uni != 0, "errorMatches", NA),
                                   # if is alt
                                   ifelse(summary_probes$count.xun != 0, "errorMatches", NA)
                                   )

# Merge
summary_probes$diagnosis<-paste(summary_probes$diagnosis1, summary_probes$diagnosis2, summary_probes$diagnosis3, sep  =",")
summary_probes$diagnosis<-gsub("NA" ,"", summary_probes$diagnosis)
summary_probes$diagnosis<-gsub(",," ,",", summary_probes$diagnosis)
summary_probes$diagnosis<-gsub("^," ,"", summary_probes$diagnosis)
summary_probes$diagnosis<-gsub(",$" ,"", summary_probes$diagnosis)

summary_probes$diagnosis1<-summary_probes$diagnosis2<-summary_probes$diagnosis3<-NULL

# Sort and clean for manual analysis
summary_probes$inversion<-substr(summary_probes$probe, 4, 12)
summary_probes<-summary_probes[order(summary_probes$inversion, summary_probes$probe), c(8, 1:7)]

# Write

write.table(summary_probes, paste0(args[3], "/readFiltering.txt"), quote = F, sep = "\t", row.names = F)

# + Summary for downloads
path<-gsub("04_qualityanalysis/","",args[3])
readcounts<-read.table(paste0(path,"01_download/readscount.txt"), sep = ",", stringsAsFactors = F)
readcounts_agg<-aggregate( V3  ~ V2, rc , sum )
colnameS(readcounts_agg)<-c("Inversion", "reads.download")
write.table(readcounts_agg, paste0(args[3], "/readcounts_aggregated.txt"), quote = F, sep = "\t", row.names = F)

#####################
# Ruth Gómez Graciani
# 07 05 2021

###############################################################################
# Description:                                                                
# Make quality control to compare known genotypes to breakseq genotypes          
###############################################################################

# LOAD ARGUMENTS 
# =========================================================================== #
args = commandArgs(trailingOnly=TRUE)

# Example
# args[1]<-"analysis/2022-11-21_1kgp_highcov_static_v2.4.1.300_v38//"  # Final read counts

# Test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("Please set input directory", call.=FALSE)
}

# CALCULATE
# =========================================================================== #

reads<-read.table(paste0(args[1],"/03_processaligned/Results_reads.txt"))

colnames(reads)<-c("inversion", "genotype", "breakpoint", "individual", "Total_reads")

reads<-reads[order(as.character(reads$genotype)),]

# Find the probability to find each genotype

      # Original strategy, all probes with the same p robability
      # fa<-read.table(paste0(args[1], "/data/datos_librerias/bplib.fa"))
      # fa<-fa[grep(">", fa$V1),]
      # fa<-as.character(fa)
      # fa<-substr(fa, 2, 13)
      # 
      # probecount<-data.frame(table(fa))
      # rownames(probecount)<-probecount$fa
      # 
      # probecount$inv<-substr(probecount$fa, 4, 13)
      # 
      # invcount<-aggregate(Freq ~ inv, probecount, sum)
      # rownames(invcount)<-invcount$inv

# New strategy, each probe has a relative weight
readfilter<-read.table(paste0(args[1],"03_processaligned/readFiltering.txt"), header = T, sep = "\t" ,fill = T)
summarized<-data.frame(fa=paste0(readfilter$orientation, readfilter$inversion), Freq=readfilter$efficiency.check, inv = readfilter$inversion)
probecount<-aggregate( Freq~ fa + inv,  summarized, sum)
rownames(probecount)<-probecount$fa

invcount<-aggregate(Freq ~ inv, probecount, sum)
rownames(invcount)<-invcount$inv

# Aggregate in general
general<-aggregate( genotype  ~ inversion + individual,  reads, paste, collapse = "/")

# # Count aggregated factors
# counted<-data.frame(table(reads$inversion, reads$individual))
# colnames(counted)<-c("inversion", "individual", "genotypes_count")

#  Take unique factors
reads_small<-reads[, c("inversion", "genotype", "individual")]
reads_small<-unique(reads_small)

# Count unique factors
counted.ag<-data.frame(table(reads_small$inversion, reads_small$individual))
colnames(counted.ag)<-c("inversion", "individual", "genotypes_aggregated_count")

# Aggregate unique
aggregated<-aggregate( genotype  ~ inversion + individual,  reads_small, paste, collapse = "/")
colnames(aggregated)<-c( "inversion",  "individual", "genotype_aggregated"  )
  
# Count reads
readscount<-aggregate( Total_reads  ~ inversion + individual,  reads, sum)


genotable<-Reduce(function(x,y) merge(x = x, y = y, all = TRUE), list(general, aggregated, counted.ag, readscount))

# Calculate error probability
  genotable$p.error<-1
  genotable[genotable$genotypes_aggregated_count > 1, "p.error"]<-0
  
  tmp<-genotable[genotable$genotypes_aggregated_count == 1 & genotable$Total_reads >1,]
  tmp$probe<-paste0(tmp$genotype_aggregated, tmp$inversion)
  
  tmp$probefreq<-probecount[tmp$probe,"Freq"]
  tmp$invcount<-invcount[tmp$inversion, "Freq"]
  
  # now the formula!
  # probability of mine  ^ probes = error calculation = probability of not detecting the other
  tmp$p.error<- (tmp$probefreq / tmp$invcount)^tmp$Total_reads
  
  genotable<-merge(genotable, tmp[,c("inversion", "individual", "p.error")], by = c("inversion", "individual"), all=TRUE)
  
  genotable$p.error<-ifelse(is.na(genotable$p.error.y), genotable$p.error.x, genotable$p.error.y)
  genotable$p.error.x<-genotable$p.error.y<-NULL
  
  # Now with het individuals
  tmp<-genotable[genotable$genotypes_aggregated_count > 1,]
  tmp$inversion<-as.character(tmp$inversion)
  tmp$individual<-as.character(tmp$individual)
  reads$individual<-as.character(reads$individual)
  reads$inversion<-as.character(reads$inversion)
  
  tmp$p.error<-apply(tmp, 1, function(x){
    # esto solo mira que haya mas de una sonda!!
    # counts<-data.frame(table(as.character(reads[reads$inversion %in%  x["inversion"] & reads$individual %in% x["individual"], "genotype"] )))
    # x["p.error"] <- ifelse(  (nrow(counts)  ==   nrow(counts[counts$Freq >= 2,])), 0, 1 )
    counts<-aggregate( Total_reads ~ genotype ,  reads[reads$inversion %in%  x["inversion"] & reads$individual %in% x["individual"],] , sum)
    x["p.error"] <- ifelse(  (nrow(counts)  ==   nrow(counts[counts$Total_reads >= 2,])), 0, 1 )
    
  })
  
  genotable<-merge(genotable, tmp[,c("inversion", "individual", "p.error")], by = c("inversion", "individual"), all=TRUE)
  genotable$p.error<-ifelse(is.na(genotable$p.error.y), genotable$p.error.x, genotable$p.error.y)
  genotable$p.error.x<-genotable$p.error.y<-NULL

  
# Save
write.table(genotable, paste0(args[1],"/03_processaligned/GTypes_FinalDataSet.txt"), col.names = T, row.names = F, quote = F, sep = "\t" )
# saveRDS(genotable, file =  paste0(args[1],"/03_processaligned/GTypes_FinalDataSet.rds"))
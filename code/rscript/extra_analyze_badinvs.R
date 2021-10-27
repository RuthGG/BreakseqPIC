rm(list=ls()) 


#####################
# Ruth Gómez Graciani
# 16 06 2021

sondas<-read.table("data/raw/datos_librerias/bplib.header.template.txt", skip = 3 , sep = "\t")

# Process & list existing probes

sondas$GENO<-substr(sondas$V2, 4, 6)
sondas$Inv<-substr(sondas$V2, 7, 15)
sondas$BP<-substr(sondas$V2, 16, 18)

####

badinvs<-read.table("data/use/invs_bad.txt", sep = "\t")
samples<-read.table("data/use/static_dataset_list.txt")
reference<-read.table("data/use/GlobalInvGenotypes_v3.1_132Invs_20210528_Genotypes.csv", sep = "\t", header = T)
genotypes<-read.table("analysis/03_processaligned/2021-05-11/min_20_30/GTypes_FinalDataSet.txt", sep = "\t", header = T)
reads<-read.table("analysis/03_processaligned/2021-05-12/min_20_30/Results_reads.txt")
reads_total<-read.table("analysis/03_processaligned/2021-06-16/min_20_0/Results_reads.txt")
colnames(reads)<-colnames(reads_total)<-c("Inv", "GENO", "BP", "Ind", "Reads")
reads$probe<-paste(reads$GENO,reads$BP , sep = "_")

reads_total$probe<-paste(reads_total$GENO,reads_total$BP , sep = "_")

# write.table(genotypes[genotypes$Inv == "HsInv0200",] , "hsinv0200_genotypes.txt", row.names = F, quote = F)

# Badinvs with bad genotypes

badinvs[badinvs$V2 == "Many incorrect genotypes",]


checkIncorrect<-function(inv){

genotypes_inv<-genotypes[genotypes$Inv == inv,c("Sample", "Inv", "Gtype_breakseq", "Total_reads")]
reference_inv<-reference[reference$Sample.ID %in% genotypes$Sample, c("Sample.ID", inv)]

comparison<-merge(genotypes_inv, reference_inv, all=T, by.x = "Sample", by.y = "Sample.ID")
comparison<-comparison[comparison$Sample %in% samples$V1,]

comparison$Result<-as.character(comparison$Gtype_breakseq) == as.character(comparison[,inv])

comparison$Change<-paste0(comparison[,inv], "to", comparison$Gtype_breakseq)

print(inv)
print(table(comparison[which(comparison$Result == "FALSE" ), "Change"]))
print(table(comparison[which(comparison$Result == "TRUE" ), "Change"]))

}

for(inv in badinvs[badinvs$V2 == "Many incorrect genotypes","V1"]){
  checkIncorrect(inv)
}

# Badinvs with few reads

badinvs[badinvs$V2 %in% c("Few reads (suspicious)", "Many individuals with 0 reads", "Many incorrect genotypes"),]

checkFew<-function(inv){

  filtered<-reads[reads$Inv == inv, ]
  total<-reads_total[reads_total$Inv == inv, ]
  

  sondas_inv<-sondas[sondas$Inv == inv, ]
  sondas_inv<-paste(sondas_inv$GENO, sondas_inv$BP , sep = "_")
  
  
  filtered$probe <- factor(filtered$probe, levels = sondas_inv)
  total$probe <- factor(total$probe, levels = sondas_inv)
  
  
 
  if(nrow(filtered)>0){
    filtered$Category<-"Filtered"
    total$Category<-"Raw"
    
      all<-rbind(filtered, total)
  }
  else{
    total$Category<-"Raw"
    all<-total
  }
  
  numbers<-data.frame(table(all$probe, all$Category))
  colnames(numbers)<-c("probe", "Category", "inds")
  
  return(
    ggplot(all)+geom_boxplot(aes(x = probe , y = Reads))+facet_grid(. ~ Category)+ggtitle(inv)+
      geom_text(data=numbers, aes(x=probe, y= Inf, vjust = 1 , label=paste0("n=",inds)))
    )
  
}

figures<-list()
for(inv in badinvs[badinvs$V2 %in% c("Few reads (suspicious)", "Many individuals with 0 reads", "Many incorrect genotypes"),"V1"]){
  figures[[inv]] <- checkFew(inv)
}

cowplot::plot_grid(plotlist =figures)


# MARIO INVS


figures<-list()
for(inv in c("HsInv1874")){
  figures[[inv]] <- checkFew(inv)
}

cowplot::plot_grid(plotlist =figures)

#  Badinvs with many reads
 

figures<-list()
for(inv in badinvs[badinvs$V2 %in% c("Many reads (suspicious)"),"V1"]){
  figures[[inv]] <- checkFew(inv)
}

cowplot::plot_grid(plotlist =figures)

# Quality check plot by the probe!

reads$inv_probe<-paste(reads$Inv, reads$probe, sep = "_")
sondas$inv_probe<-paste(sondas$Inv, sondas$GENO, sondas$BP, sep = "_")

reads$inv_probe<-factor(reads$inv_probe, levels = unique(sondas$inv_probe))

minval<-as.numeric(summary(reads$Reads)[2])

maxval<-as.numeric(summary(reads$Reads)[5])

readSummary<- aggregate(  Reads   ~ inv_probe,  reads, summary )
readSummary<-data.frame(Probe = readSummary$inv_probe, readSummary$Reads)

numbers<-data.frame(table(reads$inv_probe))
colnames(numbers)<-c("Probe", "numinds")
numbers$Probe<-factor(numbers$Probe, levels = unique(sondas$inv_probe))
readSummary<-merge(readSummary, numbers, by = "Probe", all = T)

readSummary$QualityTest<-ifelse( readSummary$numinds == 0, "No match" , 
                                 ifelse(readSummary$X3rd.Qu. < minval, "Low read count", 
                                ifelse(readSummary$X1st.Qu. > maxval , "High read count", "")))

write.csv(readSummary, "readSummary.csv")


# HsInv0200 by itself

# function to examine a single inv; needs to execute previously extra_breakseq1inv.sh
# Check that folder is correct!

OneInvAnalysis<-function(inv){

  i_reads<-read.table(paste0("analysis/02_breakseq/2021-05-11/min_20/",inv,".txt"),  sep = "\t", comment.char = "", col.names = paste0("V",seq_len(22)), fill = TRUE, quote = "", stringsAsFactors = F)
  i_ref<-read.table(paste0("analysis/02_breakseq/2021-05-11/min_20/",inv,"_refhits.txt"), sep = "\t", comment.char = "", quote = "", col.names = paste0("V",seq_len(21)), fill = T, stringsAsFactors = F)
  
  i_ref<-i_ref[which(i_ref$V3 %in% i_reads$V3),]
  refcounts<-data.frame(table(i_ref$V3))
  
  i_reads<-i_reads[,c("V1", "V2","V3", "V5")]
  i_reads<-merge(i_reads, refcounts, by.x = "V3", by.y = "Var1", all = T)
  
  i_reads$GENO<-substr(i_reads$V5, 1, 3)
  i_reads$BP<-substr(i_reads$V5, 13, 15)
  i_reads$probe<-paste(i_reads$GENO, i_reads$BP, sep = "_")
  
  # Cuantas veces mapean los reads fuera
  print("Distribution of number of matches in other parts of the genome")
  print(summary(i_reads$Freq))
  
  # Distribuión de los reads
  readsFound<-as.data.frame(table(i_reads$probe, i_reads$V1, i_reads$V2))
  print(ggplot(readsFound[readsFound$Var3 != "fil.sam",])+geom_boxplot(aes(x=Var1, y = Freq))+facet_grid(.~Var3))
  
  # Destino de los reads
  dest<-data.frame(table(i_reads$probe, i_reads$V2))
  dest<-reshape(dest, direction = "wide", idvar = "Var1", timevar = "Var2")
  dest$perc_pass<-dest$Freq.fil.sam/dest$Freq.ini.sam
  
  # Reads encontrados AFTER minqual - sin efecto de minimum length
  fil<-reads_total[reads_total$Inv == inv,]
  readsAll<-data.frame(table(fil$probe))
  colnames(readsAll)[2]<-"QualityFilter"
  dest<-merge(dest, readsAll)
  
  # Reads encontrados AFTER minqual - con efecto de minimum length
  fil<-reads[reads$Inv == inv,]
  readsUsed<-data.frame(table(fil$probe))
  colnames(readsUsed)[2]<-"ReadLengthFilter"
  dest<-merge(dest, readsUsed)
  
  dest<-dest[,c(1,3, 2, 6, 4, 5, 7, 8)]
  colnames(dest)[1]<-"Probe"
  
  return(dest)

}

OneInvAnalysis("HsInv0200")
OneInvAnalysis("HsInv1874")

# HsInv0200 with NA12004
i_reads<-read.table("analysis/02_breakseq/2021-05-11/min_20/HsInv0200.txt",  sep = "\t", comment.char = "", col.names = paste0("V",seq_len(22)), fill = TRUE, quote = "", stringsAsFactors = F)
i_ref<-read.table("analysis/02_breakseq/2021-05-11/min_20/HsInv0200_refhits.txt", sep = "\t", comment.char = "", quote = "", col.names = paste0("V",seq_len(21)), fill = T, stringsAsFactors = F)

i_ref<-i_ref[which(i_ref$V3 %in% i_reads$V3),]
refcounts<-data.frame(table(i_ref$V3))

oneind<-i_reads[i_reads$V1 =="NA12004", ]
oneind_ref<-i_ref[i_ref$V1 == "NA12004",]

write.table(oneind, "NA12004_probereads.txt", quote = F, row.names = F, col.names = F, sep = "\t")
write.table(oneind_ref, "NA12004_refreads.txt",quote = F, row.names = F, col.names = F, sep = "\t")

tb<-reads[reads$Inv == "HsInv1874",]

tb_wide<-reshape( tb[,c("Inv", "Ind", "probe","Reads" )], direction = "wide", idvar = c("Inv", "Ind"), timevar = "probe")
tb_merged<-merge(tb_wide, genotypes[genotypes$Inv == "HsInv1874",c("Inv","Sample", "STDuREF_reads", "INVuDELuINS_reads", "Gtype_breakseq")], by.x = c("Inv", "Ind"), by.y = c("Inv", "Sample"), all = T)
tb_merged[is.na(tb_merged)]<-0

write.table(tb_merged, "inv1874genotypes.txt", quote = F, sep = "\t", row.names = F)





# ----------------
sondas<-read.table("data/raw/datos_librerias/bplib.header.template.txt", skip = 3 , sep = "\t")

# Process & list existing probes

sondas$GENO<-substr(sondas$V2, 4, 6)
sondas$Inv<-substr(sondas$V2, 7, 15)
sondas$BP<-substr(sondas$V2, 16, 18)

reads<-read.table("analysis/03_processaligned/2021-05-12/min_20_30/Results_reads.txt")

colnames(reads)<-c("Inv", "GENO", "BP", "Ind", "Reads")
# reads$probe<-paste(reads$GENO,reads$BP , sep = "_")
# 
# 
# reads$inv_probe<-paste(reads$Inv, reads$probe, sep = "_")
# sondas$inv_probe<-paste(sondas$Inv, sondas$GENO, sondas$BP, sep = "_")
# 
# reads$inv_probe<-factor(reads$inv_probe, levels = unique(sondas$inv_probe))
# 
# minval<-as.numeric(summary(reads$Reads)[2])
# 
# maxval<-as.numeric(summary(reads$Reads)[5])

readSummary<- aggregate(  Reads   ~ inv_probe,  reads, summary )
readSummary<-data.frame(Probe = readSummary$inv_probe, readSummary$Reads)

numbers<-data.frame(table(reads$inv_probe))
colnames(numbers)<-c("Probe", "numinds")
numbers$Probe<-factor(numbers$Probe, levels = unique(sondas$inv_probe))
readSummary<-merge(readSummary, numbers, by = "Probe", all = T)

readSummary$QualityTest<-ifelse( readSummary$numinds == 0, "No match" , 
                                 ifelse(readSummary$X3rd.Qu. < minval, "Low read count", 
                                        ifelse(readSummary$X1st.Qu. > maxval , "High read count", "")))

write.csv(readSummary, "readSummary.csv")


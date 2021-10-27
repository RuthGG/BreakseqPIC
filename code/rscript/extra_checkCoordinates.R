

# The results from 
# run from tmp/coordCheck/
# bowtie2 -f -x ../../data/use/bowtie_index/h_sapiens_asm -U ../../data/use/bowtie_index/bplib.fa > result

samfile_original<-read.table("tmp/coordCheck/result", comment.char = "@", sep = "\t",fill = T, stringsAsFactors = F)
samfile<-samfile_original[,c(1,3,4,6)]
colnames(samfile)<-c("Probe", "Chromosome", "Coordinate", "CIGAR")

# The coordinates 

coords<-read.table("data/raw/seed_librerias/v2.3.2.300/bplib.coords", stringsAsFactors = F)
colnames(coords)<-c("Probe", "Chromosome", "Start", "End", "Length")

# Merge
check<-merge(coords, samfile, all.x = T, by = "Probe", suffixes = c("Coords", "Bowtie"))
check$Result<- ifelse(( check$ChromosomeCoords == check$ChromosomeBowtie & check$Start == check$Coordinate & check$CIGAR == "300M"), "Good", NA )

# The matches that should not be!

notREF<-samfile[ !(samfile$Probe %in% coords$Probe), ]
notREF<-notREF[!(notREF$CIGAR %in% c("*", "", " ", NA)),]

notREF<-merge(notREF, samfile_original[,c(1,5)], by.x = "Probe", by.y = "V1", all.x = TRUE)
colnames(notREF)[5]<-"MAPQ"

write.table(check, "tmp/coordCheck/REF")
write.table(notREF, "tmp/coordCheck/notREF")

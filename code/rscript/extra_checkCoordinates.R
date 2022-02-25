
###############################################################################
# Description:                                                                
# Make quality control of inversion coordinates
# After updating the library, run

# condor_submit -interactive 
# singularity shell --bind /data/bioinfo/common/bowtie2_index:/nfs/pic.es/user/r/rgomez/20210325_breakseq/data/use/bowtie_index:rw /data/bioinfo/software/rgomez_breakseq.sif
# mkdir -p 20210325_breakseq/tmp/coordCheck/
# cd  20210325_breakseq/tmp/coordCheck/
# bowtie2 -f -x ../../data/use/bowtie_index/h_sapiens_asm -U ../../data/use/bowtie_index/bplib.fa > result
# Rscript code/rscript/extra_checkCoordinates.R
###############################################################################

# Alignment
samfile_original<-read.table("coordCheck/result", comment.char = "@", sep = "\t",fill = T, stringsAsFactors = F)
samfile<-samfile_original[,c(1,3,4,6)]
colnames(samfile)<-c("Probe", "Chromosome", "Coordinate", "CIGAR")
# repeated<-samfile[samfile$Probe %in% samfile$Probe[duplicated(samfile$Probe)],]

# The coordinates
coords<-read.table(paste0("bplib.coords"), stringsAsFactors = F)
colnames(coords)<-c("Probe", "Chromosome", "Start", "End", "Length")

# Merge
check<-merge(coords, samfile, all.x = T, by = "Probe", suffixes = c("Coords", "Bowtie"))
check$Result<- ifelse(( check$ChromosomeCoords == check$ChromosomeBowtie & check$Start == check$Coordinate & check$CIGAR == "300M"), "Good", NA )

# The matches that should not be!

notREF<-samfile[ !(samfile$Probe %in% coords$Probe), ]
notREF<-notREF[!(notREF$CIGAR %in% c("*", "", " ", NA)),]

notREF<-merge(notREF, samfile_original[,c(1,5)], by.x = "Probe", by.y = "V1", all.x = TRUE)
colnames(notREF)[5]<-"MAPQ"

write.table(check, "coordCheck/REF", quote = FALSE, row.names = FALSE)
write.table(notREF, "coordCheck/notREF", quote = FALSE, row.names = FALSE)

# Interpretation

# L’arxiu REF té el resultat de comparar la llista de coordenades (pàgina bplib_hg19coords del excel de sondes) amb
# els matches en el genoma de referència d’aquelles sondes corresponents. 
# L’arxiu notREF té tots aquells matches de les sondes que no estan a la llista de coordenades (he obviat totes aquelles
# sondes que no feien match, ja que és l’esperable).

# Problemes detectables
#   
#  Sondes REF que no estan a la llista de coordenades. Això es veu perquè a la taula REF veiem les referencies en blanc o algunes no hi son
#  pero després a la notREF trobem els matches REF perduts.
# 
#  A la taula aula REF, les que tenen un match però no posa “Good” és perquè les coordenades reals no coincideixen amb les que hi ha a la taula de coordenades.
#
#  Sondes INV que fan match amb el genoma de referència, algunes millor i altres pitjor. 
#   Això vol dir que alguns reads que coincideixin amb aquestes sondes es descartaran perquè faran match en el genoma de referència, 
#   tot i que com es veu pel CIGAR i el MAPQ, seria un match molt pobre en la majoria de casos (no tots).
#   Per saber fins a quin punt això està afectant a l’anàlisi, s’hauria de comparar amb la taula del filtratge de reads.
# 

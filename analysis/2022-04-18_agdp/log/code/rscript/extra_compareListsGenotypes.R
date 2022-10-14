


gtypes<-read.table("analysis/2021-07-12_library.v2.3/03_processaligned/GTypes_FinalDataSet.txt", header= T)


gtypesInv<-read.table("analysis/2021-07-13_library.v2.3.Invs/03_processaligned/GTypes_FinalDataSet.txt",header= T)
colnames(gtypesInv)[4]<-"genotype_aggregated_inv"
gtypesOther<-read.table("analysis/2021-07-14_library.v2.3.Others/03_processaligned/GTypes_FinalDataSet.txt",header= T)
colnames(gtypesOther)[4]<-"genotype_aggregated_other"

merged<-Reduce(function(x,y) merge(x = x, y = y, by = c("inversion", "individual"), all = TRUE), 
       list(
         gtypes[c("inversion", "individual", "genotype_aggregated")],
         
         gtypesInv[c("inversion", "individual", "genotype_aggregated_inv")],
         
         gtypesOther[c("inversion", "individual", "genotype_aggregated_other")]))

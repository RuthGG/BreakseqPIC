# To know how many inds there are for each inv in the genotypes file

gfile<-read.table("data/use/GlobalInvGenotypes_v3.1_132Invs_20210528_Genotypes.csv", sep = "\t", header = T, stringsAsFactors = F)


gfile<-melt(gfile, id.vars = colnames(gfile)[1:9])

head(gfile)

gfile_complete<-gfile[!is.na(gfile$value), c("variable", "value")]

info<-data.frame(table(gfile_complete$variable))

write.table(info, file = "tmp.csv", quote = F, row.names = F, sep = "\t")

# The same but with the genotypes file from our data

gfile<-read.csv("analysis/03_processaligned/2021-05-18/smallinvs_min_20_30/GTypes_FinalDataSet.txt", header = T, sep = "\t")
# gfile<-read.table("analysis/03_processaligned/2021-04-26/min_20_30/GTypes.txt_FinalDataSet.txt", header = T)

info<-data.frame(table(gfile$Inv))
# smallinvs<-read.table("data/use/small_invs.txt")
# info<-info[as.character(info$Var1) %in% as.character(smallinvs$V1),]
write.table(info, file = "tmp.csv", quote = F, row.names = F, sep = "\t")

# ALSO, FROM A mx_plot table
# We need the function makeplot from 04_qualityanalysis_comparison.R

# Probes v2. min_20_30
args = commandArgs(trailingOnly=TRUE)

# # Example
args[1]<-"analysis/04_qualityanalysis/2021-05-18/" # Directory of results
args[2]<-"data/use/static_dataset_list.txt"  # All sample names
args[3]<-"data/use/NewInvGenotypes_v3.0_72Invs_20210426.csv" #InvFEST genotypes
condition ="smallinvs_min_20_30"
mx_plotlist<-makeplot(paste0(args[1],"/",condition,"/GTypes_compare_InvFEST.txt"),1)
mx_plot<-mx_plotlist[[1]]

mx_plot_wide<-reshape(mx_plot, direction = "wide", timevar = "Comparison", idvar = "Inversion")
# smallinvs<-read.table("data/use/small_invs.txt")
# mx_plot_wide<-mx_plot_wide[mx_plot_wide$Inversion %in% smallinvs$V1,]
write.table(mx_plot_wide, file = "tmp.csv", quote = F, row.names = F, sep = "\t")

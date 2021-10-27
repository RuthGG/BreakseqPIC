#####################
# Ruth Gómez Graciani
# 03 06 2021

###############################################################################
# Description:                                                                
# Make figures to show tagSNP analysis          
###############################################################################

# LOAD ARGUMENTS 
# =========================================================================== #
args = commandArgs(trailingOnly=TRUE)

# Example
# args[1]<-"analysis/2021-07-12_library.v2.3/04_qualityanalysis/tagSNPs_max.txt" # Input
# args[2]<-"analysis/2021-07-12_library.v2.3/04_qualityanalysis/"  # Outdir
# 

# # Test if there is at least one argument: if not, return an error
if (length(args)<2) {
  stop("One input file, one output dir", call=FALSE)
}


# LOAD PACKAGES
# =========================================================================== #

library(ggplot2)

# OPEN FILES, MAKE BASIC VARIABLES
# ======================

maxtag<-read.table(args[1], sep = "\t", header = T, stringsAsFactors = F)
originames<-colnames(maxtag)
colnames(maxtag)<-c("INV","TYPE","N_SAMPLES_F10" ,"max_SNP_10","LD_max_SNP10" )
maxtag[which(is.na(maxtag$LD_max_SNP10) | maxtag$LD_max_SNP10 %in% c("", "ND")), "LD_max_SNP10"]<-0
maxtag$LD_max_SNP10<-as.numeric(maxtag$LD_max_SNP10)

# MAKE PLOT
# ======================

plot<-ggplot(maxtag,aes(x=INV,y=LD_max_SNP10))+
  geom_point()+
  geom_hline(yintercept = 1,linetype=1)+
  geom_hline(yintercept = 0.95,linetype=2)+
  geom_hline(yintercept = 0.9,linetype=3)+
  scale_y_continuous(limits = c(0,1),breaks = c(0.7,0.8,0.9,0.95,1))+
  ylab("Linkage disequilibrium")+xlab("Inversion")+
  # theme(axis.text.x=element_text(angle = 90,  hjust = 1 ))+
  coord_flip()+
  ggtitle("Maximum LD found between each inversion and a SNP")

# SAVE FILES
# =========================================================================== #

# Save plot 
# Open file
png(paste0(args[2],"/tagSNPS_max.png"), width = 210 , height = 297 , units = "mm", res = 300)
# Create a plot
plot
# Close the file
dev.off() 
  
# MAKE GENERAL TABLE
# ======================

if (file.exists(paste0(args[2], "/genotypesCompared.txt"))){
  
  genotypesC<-read.table(paste0(args[2], "/genotypesCompared.txt"), sep ="\t", header = T)
  colnames(maxtag)<-originames
  complete<-merge(genotypesC, maxtag[,c(1,3,5)], by.x = "Inversion" , by.y = "INV", all = T)
  
  write.table(complete, file = paste0(args[2],"/genotypesCompared.txt" ), quote = F, sep = "\t", row.names = F)
  
}


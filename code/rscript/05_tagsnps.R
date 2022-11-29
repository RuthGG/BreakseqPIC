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

# # Example
# args[1]<-"analysis/2022-10-20_1kgp_highcov_static_v2.4.1.300_v38/05_tagsnps/tagSNPs_max.txt" # Input
# args[2]<-"analysis/2022-10-20_1kgp_highcov_static_v2.4.1.300_v38/05_tagsnps/"  # Outdir


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
maxtag[which(is.na(maxtag$LD_max_SNP)), ]<-0
maxtag$LD_max_SNP<-as.numeric(maxtag$LD_max_SNP)

# MAKE PLOT
# ======================

plot<-ggplot(maxtag,aes(x=INV,y=LD_max_SNP))+
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
  
rm(list=ls()) 


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


# # Example
args[1]<-"analysis/04_qualityanalysis/2021-06-15/Compare_libraries/" # Directory of results

# Test if there is at least one argument: if not, return an error
if (length(args)<1) {
  stop("One input file, one number of samples file, one output directory", call.=FALSE)
}

# LOAD PACKAGES
# =========================================================================== #
if (!require("ggplot2")) install.packages("ggplot2")
library(ggplot2)
# GENOTYPE AGREEMENT: BREAKSEQ VS. INVFEST
# =========================================================================== #


# Function that creates the table
makeplot<-function(file){
  
  # Take table
  x<-read.table(paste0(file), stringsAsFactors = FALSE, sep = "\t", header = T)
  x<-x[,c("Inversion", "Numinds.with.ref", "Numinds.geno.break", "Percentage.Different.GType", "Percentage.Filter..0.reads.", "Percentage.Same.GType", "Percentage.Filter...5.reads.")]
  
  mx_plot<-reshape2::melt(x, id.vars = c("Inversion", "Numinds.with.ref", "Numinds.geno.break"))
  names<-c("Filter (<5 reads)",
    "Filter (0 reads)",
    "Different GType",
    "Same GType")
  colors<-c("#afb6bf", "black","red","#5f9e6e")
  
  library(plyr)
  mx_plot$variable<-revalue(mx_plot$variable, c("Percentage.Different.GType"="Different GType",
                              "Percentage.Filter..0.reads."="Filter (0 reads)",
                              "Percentage.Same.GType"="Same GType",
                              "Percentage.Filter...5.reads." =  "Filter (<5 reads)"))

  # Make plot table
  mx_plot<-mx_plot[,c("Inversion","variable", "value" )]
  colnames(mx_plot)<-c("Inversion","Comparison","Genotypes compared")
  aux<-mx_plot[mx_plot$Comparison=="Same GType",]
  mx_plot$Inversion<-factor(mx_plot$Inversion,
                            levels = as.character(as.matrix(aux[order(aux$`Genotypes compared`,decreasing = T),1])))
  
  order1<-unique(as.character(as.matrix(mx_plot[mx_plot$Comparison=="Same GType",][order(as.numeric(as.matrix(mx_plot[mx_plot$Comparison=="Same GType",]$`Genotypes compared`)),decreasing = T),]$Inversion)))
  mx_plot$Inversion<-factor(mx_plot$Inversion,
                            levels = order1)
  
  mx_plot$Comparison<-factor(mx_plot$Comparison,levels = names)
  
  return(list(mx_plot, names, colors))
  
}

# Summarizing plots
condlist<-list.dirs(args[1], full.names = FALSE)
condlist<-condlist[condlist != ""]

if(length(condlist) <2){
  print("There are not multiple conditions to be compared.")
}else{
  
  mx_master<-data.frame()
  ld_master<-data.frame()
  
  for (condition in condlist) {
    mx_plotlist<-makeplot(paste0(args[1],"/",condition,"/genotypesCompared.txt"))
    
    mx_plot<-mx_plotlist[[1]]
    names<-mx_plotlist[[2]]
    colors<-mx_plotlist[[3]]
    
    mx_plot$Condition<-condition
    mx_plot<-mx_plot[!is.na(mx_plot$`Genotypes compared` ),]
    
    mx_master<-rbind(mx_master, mx_plot)
    
    ld_plot<-read.table(paste0(args[1],"/",condition,"/tagSNPs_max.txt"), sep = "\t", header = T, stringsAsFactors = F)
    # ld_plot[is.na(ld_plot)]<-0
    # ld_plot[ld_plot$LD_max_SNP10 == "ND", "LD_max_SNP10"]<-0
    ld_tab<-hist(as.numeric(ld_plot$LD_max_SNP10), breaks = c(seq(0,1,0.1)), plot = F )
    ld_tab<-data.frame(ld_tab$breaks[2:11], ld_tab$counts)
    ld_tab$Condition<-condition
    
    ld_master<-rbind(ld_master, ld_tab)
    
  }
  
  # Lineplots
  line<-ggplot(mx_master, aes(x = Inversion, y=`Genotypes compared`, group = Condition, color = Condition ))+
    geom_line()+geom_point(size = 0.5)+
    facet_grid(Comparison ~., scales = "free")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "bottom")+
    guides(color=guide_legend(ncol=3))+
    ggtitle("Per-inversion conditions comparison")
  
  # Boxplots
  box<-ggplot(mx_master)+geom_boxplot(aes(x = Condition, y =`Genotypes compared`, color = Comparison, fill = Comparison ), alpha = 0.5)+
    scale_fill_manual(values=colors,labels=names)+
    scale_color_manual(values=colors,labels=names)+
    facet_grid(Comparison ~., scales = "free")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")+
    ggtitle("Per-condition summary of results")

  boxlog<-ggplot(mx_master)+geom_boxplot(aes(x = Condition, y =log10(`Genotypes compared`), color = Comparison, fill = Comparison ), alpha = 0.5)+
    scale_fill_manual(values=colors,labels=names)+
    scale_color_manual(values=colors,labels=names)+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1), legend.position = "none")+
    facet_grid(Comparison ~., scales = "free")+
    ggtitle("Per-condition summary of results")
  
  # Tagsnps

  tags<-ggplot(ld_master)+
    geom_bar(aes(x = `ld_tab.breaks.2.11.`, y = ld_tab.counts), stat="identity")+
    geom_text(data = ld_master[ld_master$ld_tab.counts !=0,], aes(x =`ld_tab.breaks.2.11.`, y = ld_tab.counts, label =  ld_tab.counts), nudge_y = 1)+
    facet_wrap(.~ Condition)+
    ggtitle("Max LD with SNP frequency")

  
  # Open file
  png(paste0(args[1],"/comparisons.png"), width = 210 , height = 297 , units = "mm", res = 300)
  # Create a plot
  gridExtra::grid.arrange(line, box, boxlog, tags)
  # Close the file
  dev.off()
  
  
}

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

# # # # Example

# args[1]<-"analysis/2022-07-06_half1kgp_v2.3.3.300/03_processaligned/GTypes_FinalDataSet.txt" # File with results
# args[2]<-"analysis/2022-07-06_half1kgp_v2.3.3.300/data/samples.txt"  # All sample names
# args[3]<-"data/raw/GlobalInvGenotypes_v3.2_132Invs_20210528_Genotypes.csv" #InvFEST genotypes
# args[4]<- "analysis/2022-07-06_half1kgp_v2.3.3.300/04_qualityanalysis"
# args[5]<-"analysis/2022-07-06_half1kgp_v2.3.3.300/data/regions.txt"
# args[6]<-0.03 # maximum error admitted

# # Test if there is at least one argument: if not, return an error
if (length(args)<6) {
  stop("One input file, one sample names file, one reference genotypes file, one output directory, a regions list, a maximum error", call.=FALSE)
}


# LOAD PACKAGES
# =========================================================================== #

library(ggplot2)
library(reshape2)
library(grid)

# OPEN FILES, MAKE BASIC VARIABLES
# =========================================================================== #

g.ref<-read.table(args[3], sep = "\t", header = T, stringsAsFactors = F)
g.break<-read.table(args[1], header = T, stringsAsFactors = F)
inds<-read.table(args[2], stringsAsFactors = F)
invs<-read.table(args[5], stringsAsFactors = F)

# Tidy g.ref
  # Remove unused columns and rows
  inv.ref<-colnames(g.ref)[grep("Hs.*", colnames(g.ref))]
  g.ref<-g.ref[which(g.ref$Sample.ID %in% inds$V1) ,c("Sample.ID", inv.ref)]
  
  # Melt
  g.ref<-reshape2::melt(g.ref, id.vars = c("Sample.ID") )
  colnames(g.ref)<-c("Sample", "Inv", "Gtype_experimental")
  
  # Filter g.ref
  g.ref[which(g.ref$Gtype_experimental == "ND"), "Gtype_experimental"]<-NA
  g.ref$Inv <- as.character(g.ref$Inv)
  g.ref<-g.ref[which(g.ref$Inv %in% invs$V1),]
  
  # Make all caps
  g.ref$Gtype_experimental<-toupper(g.ref$Gtype_experimental)

# Filter g.break
  colnames(g.break)[c(1,2)]<-c("Inv", "Sample")
  g.break<-g.break[which(g.break$Inv %in% invs$V1),]
  g.break.original<-g.break

# HOW MANY INDIVIDUALS ARE GENOTYPED FOR EACH INV
# =========================================================================== #

# In experimental data
complete.cases<-g.ref[!is.na(g.ref$Gtype_experimental), ]
numinds.ref<-data.frame(table(complete.cases$Inv))
colnames(numinds.ref)<-c("Inv", "Numinds.geno.ref")

# In Breakseq data
numinds.break<-data.frame(table(g.break$Inv))
colnames(numinds.break)<-c("Inv", "Numinds.geno.break")


# GENOTYPE REARRANGEMENT TO MAKE THEM COMPARABLE
# =========================================================================== #

# Sort  genotypes alphabetically (breakseq is already alphabetical)
# Experimental
g.ref.split<-cbind(g.ref[,c("Sample", "Inv")], data.frame(do.call("rbind", strsplit(as.character(g.ref$Gtype_experimental), "/", fixed = TRUE))))
g.ref.split$X1 <- as.character(g.ref.split$X1 )
g.ref.split$X2 <- as.character(g.ref.split$X2 )

g.ref.split<-melt(g.ref.split, id.vars = c("Sample", "Inv"), )
g.ref.split<-g.ref.split[order(g.ref.split$value), ]

g.ref<-aggregate( value  ~ Inv + Sample,  g.ref.split, paste, collapse = "/")
colnames(g.ref)<-c("Inv", "Sample", "Gtype_experimental")

counted<-data.frame(table(g.ref.split$Inv, g.ref.split$Sample))
colnames(counted)<-c("Inv", "Sample","genotypes_ref_count")

g.ref<-merge(g.ref, counted)


# GENOTYPE AGREEMENT: BREAKSEQ VS. INVFEST
# =========================================================================== #

# Make comparison
comparison<-merge(g.ref, g.break, all=T)

# Expand genotypes when necessary
# To make them comparable, even X males are made diploid
comparison$pastedref<-paste(comparison$Gtype_experimental, comparison$Gtype_experimental, sep = "/")
comparison[which(comparison$genotypes_ref_count == 1),"Gtype_experimental"]<-comparison[which(comparison$genotypes_ref_count == 1),"pastedref"]

comparison$pasted<-paste(comparison$genotype_aggregated, comparison$genotype_aggregated, sep = "/")
comparison$Gtype_breakseq <- comparison$genotype_aggregated
comparison[which(comparison$genotypes_aggregated_count == 1),"Gtype_breakseq"]<-comparison[which(comparison$genotypes_aggregated_count == 1),"pasted"]

# Classify cases
comparison$class<- ifelse(is.na(comparison$Gtype_experimental), "No Reference",
                          ifelse(is.na(comparison$Total_reads),  "Filter (0 reads)", 
                                 ifelse(comparison$p.error <= args[6] , 
                                        ifelse(comparison$Gtype_experimental == comparison$Gtype_breakseq,
                                               "Same GType", "Different GType"),
                                        "Filter (few reads)" )))

# PLOT OF THE COMPARISON
# =========================================================================== #

comp.withref<-comparison[comparison$class != "No Reference", ]
counts<-data.frame(table(comp.withref$Inv, comp.withref$class))
colnames(counts)<-c("Inv", "class", "count")
percentage<-merge(counts, numinds.ref)
percentage$percentage<-percentage$count / percentage$Numinds.geno.ref*100

names<-c("Filter (0 reads)",
         "Filter (few reads)",
         "Different GType",
         "Same GType")
colors<-c("black","#afb6bf", "red","#5f9e6e")

# Sort factors
colnames(percentage)<-c("Inversion","Comparison","Count", "Numinds.geno.ref", "Genotypes compared")
aux<-percentage[percentage$Comparison=="Same GType",]
percentage$Inversion<-factor(percentage$Inversion,
                          levels = as.character(as.matrix(aux[order(aux$`Genotypes compared`,decreasing = T),1])))

order1<-unique(as.character(as.matrix(percentage[percentage$Comparison=="Same GType",][order(as.numeric(as.matrix(percentage[percentage$Comparison=="Same GType",]$`Genotypes compared`)),decreasing = T),]$Inversion)))
percentage$Inversion<-factor(percentage$Inversion,
                          levels = order1)

percentage$Comparison<-factor(percentage$Comparison,levels = names)

# Take only necessary names and colors
colors<-colors[names %in% unique(percentage$Comparison)]
names<-names[names %in% unique(percentage$Comparison)]

# Plot 1
plot_compare<-ggplot(percentage) +
  geom_bar(aes(Inversion,
               as.numeric(as.matrix(`Genotypes compared`)),
               fill=Comparison),position="stack", stat="identity") +
  scale_fill_manual(values=colors,labels=names)+
  xlab("Inversion")+ylab("Genotype agreement")+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "top")+
  ggtitle("Comparison of Breakseq results with known genotypes")

# Plot 2
  comp.withref$Prediction<-ifelse(comp.withref$Gtype_breakseq == comp.withref$Gtype_experimental, "Correct", "Different")
  comp.withref<-comp.withref[!is.na(comp.withref$Prediction), ]
  plot_readDistribution<-ggplot(comp.withref, aes(x = Prediction, y = p.error))+
    geom_boxplot(outlier.alpha = 0)+
    geom_jitter(width = 0.25, alpha = 0.15)+
    ggtitle("Accuracy compared to distribution of error probabilities")
  
# Plot 2.1
  readfilter<-c(0, 0.001, 0.01, 0.02, 0.03, 0.04, 0.05)
  goodbad<-data.frame()
  comp.withref$Prediction<-as.factor(comp.withref$Prediction)
  for (n in readfilter) {
    
    name<-n
    
    tmp<-comp.withref[which(comp.withref$p.error <= name),]
    tmp_bis<-data.frame(table(tmp$Prediction))
    goodbad<-rbind(goodbad, data.frame(f=paste0(name) , tmp_bis))
  }
  
  goodbad<-reshape(goodbad,idvar="f", timevar = "Var1", direction = "wide")  
  goodbad$total<-goodbad$Freq.Correct + goodbad$Freq.Different

  goodbad<-melt(goodbad, id.vars = c("f", "total"), variable.name = "Count", value.name = "Count.value")
  goodbad$Count.percentage<- goodbad$Count.value/goodbad$total
  
  filtercheck<-ggplot(goodbad)+geom_bar( aes(x = as.factor(f),  y  = Count.value, fill = Count ), stat = "identity", position = "dodge")+
    geom_text ( aes(x = as.factor(f), label= paste0(Count.value, "\n (", round(Count.percentage*100, 2), "%)"), y  = 1, group = Count ),
                position = position_dodge(width = 1), vjust = -0.1)+
    facet_grid(Count  ~.,scales = "free")+
    scale_fill_manual(values = c("#5f9e6e","red"))+
    theme(legend.position = "none")+
    ggtitle("Percentage of successes and errors depending on filter")
  
# Plot 3
  g.b.summary<-aggregate( Total_reads  ~ Inv,g.break.original, summary)

  g.b.summary<-data.frame(Inv = g.b.summary$Inv, g.b.summary$Total_reads)
 
  g.break.original$Inv<-factor(g.break.original$Inv , levels=  g.b.summary$Inv[order(g.b.summary$Max.)] ) 
  
  
  minval<-as.numeric(summary(g.break.original$Total_reads)[2])
  
  maxval<-as.numeric(summary(g.break.original$Total_reads)[5])
  

  readDist<-ggplot(g.break.original)+geom_boxplot(aes(x = Inv, y = Total_reads))+coord_flip()+
    annotate("rect", xmin=-Inf, xmax=Inf, ymin = minval , ymax = maxval, alpha=0.2, fill="yellow")+
     ggtitle("Per-inversion read distribution")

# Plot with changes
  cs<-comparison[which(!is.na(comparison$Gtype_experimental) & !is.na(comparison$Gtype_breakseq)) ,c("Inv", "Sample", "Gtype_experimental", "Gtype_breakseq", "p.error", "class")]
  cs<-cs[cs$class == "Different GType",]
  cs$change<-paste0(cs$Gtype_experimental, "->", cs$Gtype_breakseq)  
  changes<-data.frame(table(cs$Inv, cs$change)  )
  
  
  cplot<-ggplot(changes[changes$Freq > 0,],aes(x = Var1, y = Freq, fill=Var2))+geom_bar(, stat = "identity")+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
    geom_text(aes(label=Freq), position = position_stack(vjust= 0.5),
              colour = "white")+
    ggtitle("Direction of changes for each inversion")
  

# SAVE FILES
# =========================================================================== #

# Save plot 
  # Open file
  png(paste0(args[4],"/genotypesCompared.png"), width = 210 , height = 297 , units = "mm", res = 300)
  # Create a plot
  # ggpubr::ggarrange( plot_compare, plot_readDistribution)
  grid.newpage()
  # Create layout : nrow = 3, ncol = 2
  pushViewport(viewport(layout = grid.layout(nrow = 3, ncol = 1)))
  # A helper function to define a region on the layout
  define_region <- function(row, col){
    viewport(layout.pos.row = row, layout.pos.col = col)
  } 
  # Arrange the plots
  print(plot_compare, vp = define_region(row = 1, col = 1))   
  print(cplot, vp = define_region(row = 2, col = 1))
  print(filtercheck, vp = define_region(row = 3, col = 1))
  # Close the file
  dev.off() 
  
  # Open file
  png(paste0(args[4],"/generalReadDistribution.png"), width = 210 , height = 297 , units = "mm", res = 300)
  # Create a plot
  readDist
  
  # Close the file
  dev.off() 

# Save data
  colnames(percentage)<-c("Inversion", "Category", "Count", "Numinds.with.ref", "Percentage")
  
  percentage.wide<- reshape(percentage, drop = "Count",v.names= c("Percentage"), timevar = "Category", idvar = c("Inversion", "Numinds.with.ref") , direction = "wide")
  percentage.wide.complete<-merge(percentage.wide, numinds.break, by.x = "Inversion", by.y = "Inv", all =TRUE)
  # This last one is just in case...
  percentage.wide.complete<-merge(percentage.wide.complete, invs, by.x = "Inversion", by.y = "V1", all=T)
  
  # Sort
  percentage.wide.complete$Inversion <- as.character(percentage.wide.complete$Inversion )
  percentage.wide.complete<- percentage.wide.complete[order(percentage.wide.complete$Inversion ),]
  percentage.wide.complete[is.na(percentage.wide.complete$Numinds.geno.break),"Numinds.geno.break" ]<-0
  
  # Add info about reads
  g.b.summary$QualityTest<-ifelse(g.b.summary$X3rd.Qu. < minval, "Low read count", 
                                                 ifelse(g.b.summary$X1st.Qu. > maxval , "High read count", ""))
  
  
  percentage.wide.complete<-merge(percentage.wide.complete, g.b.summary[,c("Inv", "Min.", "Median", "Mean", "Max.", "QualityTest")], by.x = "Inversion", by.y = "Inv")
  colnames(percentage.wide.complete)[colnames(percentage.wide.complete) %in% c("Min.", "Median", "Mean", "Max.")]<-c("Reads.Min", "Reads.Median", "Reads.Mean", "Reads.Max")
  

  write.table(percentage.wide.complete, file = paste0(args[4],"/genotypesQualityReport.txt" ), quote = F, sep = "\t", row.names = F)

  write.table(comparison[,c("Inv", "Sample", "Gtype_breakseq", "genotypes_aggregated_count","Total_reads", "p.error")],  file = paste0(args[4],"/genotypesClean.txt" ), quote = F, sep = "\t", row.names = F)

  write.table(comparison[which(!is.na(comparison$Gtype_experimental) & !is.na(comparison$Gtype_breakseq)) ,c("Inv", "Sample", "Gtype_experimental", "Gtype_breakseq", "p.error", "class")], file = paste0(args[4],"/genotypesCompared.txt" ), quote = F, sep = "\t", row.names = F) 
  
  write.table(changes[changes$Freq >0,], file = paste0(args[4],"/genotypesComparedChanges.txt" ), quote = F, sep = "\t", row.names = F) 
  
 
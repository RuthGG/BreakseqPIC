rm(list=ls()) 


############
## GTYPES ##
############
x<-read.table("/home/jon/Desktop/Breakseq_allINVs/Discussion/Analysis_InvFEST_genotypes/GTypes_compare_InvFEST.txt")

KGP<-read.table("/home/jon/Desktop/Functional_Analysis/GEUVADIS/datos/1000GP_info")
invfest=read.table("/home/jon/Desktop/Breakseq_allINVs/Data/InvFEST_GTypes/InvGenotypes_v4.8_45invs_complete_20170510_UnrelatedHapMap.txt",
                   sep="\t")
common<-intersect(as.character(as.matrix(invfest$V1)),as.character(as.matrix(KGP$V1)))
common<-length(common)

invs<-unique(as.character(as.matrix(x$V2)))
column1<-c()
column2<-c()
column3<-c()
for(inv in invs){
  
  x_esp<-x[x$V2==inv,]
  v1<-as.character(as.matrix(x_esp[x_esp$V3>4,6]))
  v2<-as.character(as.matrix(x_esp[x_esp$V3>4,7]))
  
  bien<-as.numeric(table(v1==v2)["TRUE"])
  mal<-as.numeric(table(v1==v2)["FALSE"])
  filtered<-dim(x_esp)[1]-dim(x_esp[x_esp$V3>4,])[1]
  no_reads<-common-dim(x_esp)[1]
  
  column1<-c(column1,rep(inv,4))
  column2<-c(column2,c("Same GType","Different GType","Filter (<5 reads)","Filter (0 reads)"))
  column3<-c(column3,c(bien,mal,filtered,no_reads))
}

df<-data.frame(inv=column1,
               comparison=column2,
               numbers=column3)

library(ggplot2)
mx_plot<-df

colnames(mx_plot)<-c("Inversion","Comparison","Genotypes compared")
aux<-mx_plot[mx_plot$Comparison=="Same GType",]
mx_plot$Inversion<-factor(mx_plot$Inversion,
                          levels = as.character(as.matrix(aux[order(aux$`Genotypes compared`,decreasing = T),1])))

order1<-unique(as.character(as.matrix(mx_plot[mx_plot$Comparison=="Same GType",][order(as.numeric(as.matrix(mx_plot[mx_plot$Comparison=="Same GType",]$`Genotypes compared`)),decreasing = T),]$Inversion)))
mx_plot$Inversion<-factor(mx_plot$Inversion,
                          levels = order1)

mx_plot$Comparison<-factor(mx_plot$Comparison,levels = c("Filter (0 reads)",
                                                         "Filter (<5 reads)",
                                                         "Different GType",
                                                         "Same GType"))

ggplot(mx_plot) +
  geom_bar(aes(Inversion,
               as.numeric(as.matrix(`Genotypes compared`)),
               fill=Comparison),position="stack", stat="identity") +
  scale_fill_manual(values=c("black","#afb6bf", "red","#5f9e6e"),
                    labels=c("Filter (0 reads)",
                             "Filter (<5 reads)",
                             "Different GType",
                             "Same GType"))+
  xlab("Inversion")+ylab("Genotype agreement")+
  scale_y_continuous(limits = c(0,434),breaks = c(0,200,400))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "top")


########
## LD ##
########
dir<-"/home/jon/Desktop/Breakseq_allINVs/Discussion/Linkage_disequilibrium_new/"
ld<-read.table("/home/jon/Desktop/Functional_Analysis/Extended_Functional/GWAS/Tabla_LD_all/GLOBAL_TagSNPs_INVs_v1.0_GlobalPaper",
               header=T)

column1<-c()
column2<-c()
column3<-c()
for(inv in invs){
  y<-read.table(paste0(dir,"/",inv,"/plink_",inv,"_filter5.ld"),header=T)
  y<-y[y$SNP_B!=inv,]
  y<-y[grep("rs",y$SNP_B),]
  y_esp<-head(y[order(y$R2,decreasing = T),],1)
  
  ld_esp<-ld[ld$Inv==inv,]
  ld_snps<-as.character(as.matrix(ld_esp[ld_esp$GLB==1,4]))
  
  if(length(na.omit(match(as.character(as.matrix(y_esp$SNP_B)),ld_snps)))!=0){
    print(inv)
    print("coincide")
    column1<-c(column1,inv)
    column2<-c(column2,as.character(as.matrix(y_esp$SNP_B)))
    column3<-c(column3,as.numeric(as.matrix(y_esp$R2)))
    
  }else{
    ld_esp<-ld_esp[grep("rs",ld_esp$SNP),]
    ld_snps<-as.character(as.matrix(head(ld_esp[order(ld_esp$GLB,ld_esp$EUR,decreasing = T),],1)[,4]))
    if(ld_snps==(as.character(as.matrix(y_esp$SNP_B)))){
      print(inv)
      print("coincide")
      column1<-c(column1,inv)
      column2<-c(column2,as.character(as.matrix(y_esp$SNP_B)))
      column3<-c(column3,as.numeric(as.matrix(y_esp$R2)))
      
    }else{
      if(inv=="HsInv0102"){
        print(inv)
        print("coincide")
        column1<-c(column1,inv)
        column2<-c(column2,as.character(as.matrix(y_esp$SNP_B)))
        column3<-c(column3,as.numeric(as.matrix(y_esp$R2)))
        
      }else{
        print(inv)
        print("FAIL")
        break
      }

    }
  }
}
df<-data.frame(inv=column1,
               tagsnp=column2,
               r2=column3)


df$inv<-factor(df$inv,
               levels=as.character(as.matrix(aux[order(aux$`Genotypes compared`,decreasing = T),1])))
mx_plot$Inversion<-factor(mx_plot$Inversion,
                          levels = as.character(as.matrix(aux[order(aux$`Genotypes compared`,decreasing = T),1])))


ggplot(df,aes(x=inv,y=r2))+
  geom_point(size=3)+
  geom_hline(yintercept = 1,linetype=1)+
  geom_hline(yintercept = 0.95,linetype=2)+
  geom_hline(yintercept = 0.9,linetype=3)+
  theme_classic()+
  scale_y_continuous(limits = c(0.6,1),breaks = c(0.7,0.8,0.9,0.95,1))+
  ylab("Linkage disequilibrium")
  

#################
## MULTIPROBES ##
#################
inv_1306<-read.table("/home/jon/Desktop/Breakseq_allINVs/Discussion/Sondas_multiples/GTypes_sondas_multiples_ALT_HsInv1306.txt")
inv_1141<-read.table("/home/jon/Desktop/Breakseq_allINVs/Discussion/Sondas_multiples/GTypes_sondas_multiples_HsInv1141.txt")


# inv_1306
inv_1306<-inv_1306[inv_1306$V3>4 & inv_1306$V7>4,]
inv_1306<-inv_1306[inv_1306$V6!="ND",]
inv_1306<-inv_1306[inv_1306$V10!="ND",]

inv_1306<-inv_1306[,c(6,10)]
tt=table(as.character(as.matrix(inv_1306$V6)),as.character(as.matrix(inv_1306$V10)))


library(reshape2)
co=melt(tt)

ggplot(co, aes(Var1, Var2)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = value))+
  geom_text(aes(fill = co$value, label = round(co$value, 2))) + # write the values
  scale_fill_gradient2(low = ("darkred"), 
                       mid = "white", 
                       high = ("midnightblue"), 
                       midpoint = 0)+
  theme(panel.grid.major.x=element_blank(), #no gridlines
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"), # background=white
        axis.text.x = element_text(angle=90, hjust = 1,vjust=1,size = 12,face = "bold"),
        plot.title = element_text(size=20,face="bold"),
        axis.text.y = element_text(size = 12,face = "bold"))+
  scale_x_discrete(name="Set1") +
  scale_y_discrete(name="Set2") +
  ggtitle("HsInv1306") + 
  theme(legend.title=element_text(face="bold", size=14)) +
  labs(fill="Number of GTypes")


# inv_1141
inv_1141<-inv_1141[inv_1141$V3>4 & inv_1141$V7>4,]
inv_1141<-inv_1141[inv_1141$V6!="ND",]
inv_1141<-inv_1141[inv_1141$V10!="ND",]

inv_1141<-inv_1141[,c(6,10)]
tt=table(as.character(as.matrix(inv_1141$V6)),as.character(as.matrix(inv_1141$V10)))


library(reshape2)
co=melt(tt)

ggplot(co, aes(Var1, Var2)) + # x and y axes => Var1 and Var2
  geom_tile(aes(fill = value))+
  geom_text(aes(fill = co$value, label = round(co$value, 2))) + # write the values
  scale_fill_gradient2(low = ("darkred"), 
                       mid = "white", 
                       high = ("midnightblue"), 
                       midpoint = 0)+
  theme(panel.grid.major.x=element_blank(), #no gridlines
        panel.grid.minor.x=element_blank(), 
        panel.grid.major.y=element_blank(), 
        panel.grid.minor.y=element_blank(),
        panel.background=element_rect(fill="white"), # background=white
        axis.text.x = element_text(angle=90, hjust = 1,vjust=1,size = 12,face = "bold"),
        plot.title = element_text(size=20,face="bold"),
        axis.text.y = element_text(size = 12,face = "bold"))+
  scale_x_discrete(name="Set1") +
  scale_y_discrete(name="Set2") +
  ggtitle("HsInv1141") + 
  theme(legend.title=element_text(face="bold", size=14)) +
  labs(fill="Number of GTypes")












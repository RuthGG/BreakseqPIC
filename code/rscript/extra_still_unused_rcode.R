
################################################################
## 1) LINKAGE DISEQUILIBRIUM COMPARISON WITH 1KGP ##
################################################################
  
  ##### LD with BREAKSEQ
  # dir<-"/home/jon/Desktop/Breakseq_allINVs/Discussion/Linkage_disequilibrium_new/"
  # column1<-c()
  # column2<-c()
  # column3<-c()
  # for(ff in folders){
  #     
  #   if(file.exists(paste0(dir,ff,"/plink_",ff,"_filter5.ld"))){
  #     x<-read.table(paste0(dir,ff,"/plink_",ff,"_filter5.ld"),header=T)
  #     x<-x[x$SNP_B!=ff,]
#     ld<-max(x[order(x$R2,decreasing = T),]$R2)
#     
#     column1<-c(column1,ff)
#     column2<-c(column2,ld)
#   }
#   # if(file.exists(paste0(dir,ff,"/diagnostico_5.txt"))){
#   #   x<-read.table(paste0(dir,ff,"/diagnostico_5.txt"),header=T)
#   #   x<-x[x$SNP_B!=ff,]
#   #   ld<-max(x[order(x$R2,decreasing = T),]$R2)
#   #   
#   #   column1<-c(column1,ff)
#   #   column2<-c(column2,ld)
#   # }
#   
# 
# }
# df<-data.frame(inv=column1,ld_bs=column2)
### BREAKSEQ DATA
# dir<-"analysis/05_tagsnps/2021-05-12/min_20_30/"
# dir<-"analysis/05_tagsnps/2021-05-04/all/"
# breakseq<-read.table(paste0(dir, "summary_tagSNPs.txt"), header = TRUE)
# 
# df<-breakseq[c("INV", "LD_10")]
# colnames(df)<-c("inv", "ld_bs")
# df<-df[df$ld_bs != "<0.5",]
# df$ld_bs <- as.numeric(as.character(df$ld_bs))
# 
# plot<-ggplot(df)+geom_point(aes(x = inv, y = ld_bs))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
#   geom_hline(yintercept=0.9, color = "blue" )+geom_hline(yintercept=0.95, color = "red" )+geom_hline(yintercept=1 )
# 
# 
# # Open file -- Está apaisado!!
# png(paste0(dir, "linkage_disequilibrium.png") , height = 210 , width = 297 , units = "mm", res = 300)
# # Create a plot
# plot
# # Close the file
# dev.off() 
# 
# ##### LD with 1KGP
# dir<-"../data/use/1KGPtagSNPs_small/"
# folders<-list.dirs(dir, full.names=FALSE)
# folders<-folders[grep("HsInv",folders)]
# column1<-c()
# column2<-c()
# column3<-c()
# for(ff in folders){
#   print(ff)
#   if(file.exists(paste0(dir,ff,"/plink.ld"))){
#     x<-read.table(paste0(dir,ff,"/plink.ld"),header=T, stringsAsFactors = FALSE)
#     inv<-unique(as.character(as.matrix(x$SNP_A)))
#     x<-x[x$SNP_B!=inv,]
#     # x<-x[grep("rs",x$SNP_B),]
#     ld<-max(x[order(x$R2,decreasing = T),]$R2)
#     
#     
#     column1<-c(column1,ff)
#     column2<-c(column2,ld)
#   }else if(length(dir(paste0(dir,ff)))>1){
#     coords=read.table("../data/raw/datos_librerias/bplib.coords")
#     if(ff=="HsInv0052_region"){
#       coords<-c(coords[coords$V1=="REFHsInv0052BP1HG18",3],coords[coords$V1=="REFHsInv0052BP1HG18",4])
#     }else{
#       coords<-coords[grep(ff,coords$V1),][1,c(3,4)]
#     }
#     
#     
#     counter<-1
#     for(file in dir(paste0(dir,ff))){
#       if(counter==1){
#         dfx<-read.table(paste0(dir,ff,"/",file),header=T)
#         counter<-counter+1
#       }else{
#         x<-read.table(paste0(dir,ff,"/",file),header=T)
#         dfx<-rbind(dfx,x)
#         counter<-counter+1
#       }
#     }
#     
#     select_dfx<-unique(dfx[,c(2,3)])
#     select_dfx<-as.character(as.matrix(select_dfx[select_dfx$BP_A>(as.numeric(coords[1])-10000) & select_dfx$BP_A<(as.numeric(coords[2])+10000),2]))
#     if(ff=="HsInv0052_region"){
#       select_dfx<-"BI_GS_CNV_3_162512202_162525833"
#     }
#     dfx<-dfx[dfx$SNP_A==select_dfx,]
#     
#     x<-dfx
#     inv<-unique(as.character(as.matrix(x$SNP_A)))
#     x<-x[x$SNP_B!=inv,]
#     x<-x[grep("rs",x$SNP_B),]
#     ld<-max(x[order(x$R2,decreasing = T),]$R2)
#     
#     
#     column1<-c(column1,ff)
#     column2<-c(column2,ld)
#   }
#   
#   
# }
# df2<-data.frame(inv=column1,ld_kgp=column2)
# df2<-df2[df2$ld_kgp > 0,]
# df2<-aggregate(  ld_kgp ~ inv , df2, max)
# 
# ##### PLOT COMPARISON
# not_detected<-setdiff(df$inv,df2$inv)
# common<-intersect(df$inv,df2$inv)
# rownames(df)<-df$inv
# rownames(df2)<-df2$inv
# df_c<-df[common,]
# df2_c<-df2[common,]
# df_c$ld_kgp<-df2_c$ld_kgp
# 
# df_c$ld_bs<-1-df_c$ld_bs
# df_c$ld_kgp<-1-df_c$ld_kgp
# 
# df_plot<-data.frame(inv=c(as.character(as.matrix(df_c$inv)),as.character(as.matrix(df_c$inv))),
#                     Approach=c(rep("BS",dim(df_c)[1]),rep("KGP",dim(df_c)[1])),
#                     ld=c(df_c$ld_bs,df_c$ld_kgp))
# df_plot$inv<-as.character(as.matrix(df_plot$inv))
# df_plot$Approach<-as.character(as.matrix(df_plot$Approach))
# df_plot$ld<-as.numeric(as.matrix(df_plot$ld))
# 
# df_plot<-rbind(df_plot,c("Mean","BS_mean",mean(df_c$ld_bs)))
# df_plot<-rbind(df_plot,c("Mean","KGP_mean",mean(df_c$ld_kgp)))
# df_plot$Approach<-as.factor(df_plot$Approach)
# df_plot$ld<-as.numeric(as.matrix(df_plot$ld))
# 
# df_plot_points<-data.frame(inv=df_c$inv,diff=c(df_c$ld_bs-df_c$ld_kgp),Approach="BS")
# order<-as.character(as.matrix(df_plot_points[order(df_plot_points$diff),1]))
# df_plot$inv<-factor(df_plot$inv,levels=c(order,"Mean"))
# 
# ggplot(df_plot, aes(x = as.factor(inv), 
#                         y = ld * ((-1)^(Approach == "KGP" | Approach == "KGP_mean")), 
#                         fill = Approach)) + 
#   geom_bar(stat = "identity") +
#   scale_x_discrete(name = "Inversion") +
#   #geom_point(data=df_plot_points,aes(as.factor(inv),diff),col="black")+
#   scale_y_continuous(name = "1-LD", breaks=c(-0.4, -0.2, 0,0.2,0.4),
#                      limits = c(-0.5,0.5))+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         legend.position = "top")+
#   geom_hline(yintercept = 0,linetype=3)+
#   scale_fill_manual(values=c("chocolate1","chocolate3", "steelblue2","steelblue4"))
#                     
# 
# dim(df_plot_points[df_plot_points$diff>0,])
# dim(df_plot_points[df_plot_points$diff<0,])
# dim(df_plot_points[df_plot_points$diff==0,])
# 
# 
# #### Plot not detected
# df_plot2<-df[not_detected,]
# df_plot2$Approach<-"BS"
# colnames(df_plot2)<-c("inv","ld","Approach")
# df_plot2<-df_plot2[,c(1,3,2)]
# df_plot2$ld<-1-df_plot2$ld
# 
# df_plot2$inv<-as.character(as.matrix(df_plot2$inv))
# df_plot2$Approach<-as.character(as.matrix(df_plot2$Approach))
# df_plot2$ld<-as.numeric(as.matrix(df_plot2$ld))
# df_plot2<-rbind(df_plot2,c("Mean","BS_mean",mean(df_plot2$ld)))
# df_plot2$Approach<-as.factor(df_plot2$Approach)
# df_plot2$ld<-as.numeric(as.matrix(df_plot2$ld))
# 
# ggplot(df_plot2, aes(x = as.factor(inv), 
#                     y = ld, 
#                     fill = Approach)) + 
#   geom_bar(stat = "identity") +
#   scale_x_discrete(name = "Inversion") +
#   scale_y_continuous(name = "1-LD", breaks=c(-0.4, -0.2, 0,0.2,0.4),
#                      limits = c(-0.5,0.5))+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         legend.position = "top")+
#   geom_hline(yintercept = 0,linetype=3)+
#   scale_fill_manual(values=c("chocolate1","chocolate3", "steelblue2","steelblue4"))+
#   annotate("text", x = 7, y = -0.2, label = "Not detected",size=7)
# 
# 
# ## Boxplot
# df_plot3<-df_plot[df_plot$Approach=="BS" | df_plot$Approach=="KGP",]
# ggplot(df_plot3, aes(y = Approach, 
#                     x = ld,fill=Approach))+ 
#   geom_boxplot()+
#   scale_x_continuous(name = "1-LD")+
#   theme_classic()+
#   scale_fill_manual(values=c("chocolate1", "steelblue2"))
#   
# options(warn=1)
# wilcox.test(df_plot3[df_plot3$Approach=="BS",]$ld,df_plot3[df_plot3$Approach=="KGP",]$ld,
#             exact = FALSE,alternative = "two.sided")
# 
# 
# ## Join plot1 and 2
# df_plot4<-rbind(df_plot,df_plot2)
# df_plot4$Group<-c(rep("Detected",dim(df_plot)[1]),rep("Not detected",dim(df_plot2)[1]))
# 
# ggplot(df_plot4, aes(x = (inv), 
#                     y = ld * ((-1)^(Approach == "KGP" | Approach == "KGP_mean")), 
#                     fill = Approach)) + 
#   geom_bar(stat = "identity") +
#   scale_x_discrete(name = "Inversion") +
#   #geom_point(data=df_plot_points,aes(as.factor(inv),diff),col="black")+
#   scale_y_continuous(name = "1-LD", breaks=c(-0.4, -0.2, 0,0.2,0.4),
#                      limits = c(-0.5,0.5))+
#   theme_classic()+
#   theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
#         legend.position = "top")+
#   geom_hline(yintercept = 0,linetype=3)+
#   scale_fill_manual(values=c("chocolate1","chocolate3", "steelblue2","steelblue4"))+
#   facet_grid(. ~ Group, scales = "free_x", space = "free_x")




################################################################
## 2) LINKAGE DISEQUILIBRIUM WITH KNOWN TAG SNPS FROM INVFEST ##
################################################################
## Miramos si el tag SNP reportado con el BreakSeq es el mismo que con los datos experimentales del NatComm
dir<-paste0("/home/jon/Desktop/BreakSeq_FunctionalPaper/Analysis/Discussion/",ana_f,"/Linkage_disequilibrium/")
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
  



########################################
## 3) INVERSIONS WITH MULTIPLE PROBES ##
########################################
dir<-paste0("/home/jon/Desktop/BreakSeq_FunctionalPaper/Analysis/Discussion/",ana_f,"/Multiple_probes/")

inv_1306<-read.table(paste0(dir,"/GTypes_sondas_multiples_ALT_HsInv1306.txt"))
inv_1141<-read.table(paste0(dir,"/GTypes_sondas_multiples_HsInv1141.txt"))


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








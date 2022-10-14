 
rm(list=ls())
options(warn=2)
library(ggplot2)
# folders<-dir("/home/jon/Desktop/Breakseq_allINVs/Discussion/Linkage_disequilibrium_new/")
# folders<-folders[grep("HsInv",folders)]
# exclude<-c("HsInv1275","HsInv1819","HsInv1794","HsInv0901","HsInv1848","HsInv1468","HsInv1070","HsInv0102")
# folders<-setdiff(folders,exclude)


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

breakseq<-read.table("../analysis/05_tagsnps/2021-05-04/all/summary_tagSNPs.txt", header = TRUE)

df<-breakseq[c("INV", "LD_10")]
colnames(df)<-c("inv", "ld_bs")
df<-df[df$ld_bs != "<0.5",]
df$ld_bs <- as.numeric(as.character(df$ld_bs))

ggplot(df)+geom_point(aes(x = inv, y = ld_bs))+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_hline(yintercept=0.9, color = "blue" )+geom_hline(yintercept=0.95, color = "red" )+geom_hline(yintercept=1 )

##### LD with 1KGP
dir<-"../data/use/1KGPtagSNPs_small/"
folders<-list.dirs(dir, full.names=FALSE)
folders<-folders[grep("HsInv",folders)]
column1<-c()
column2<-c()
column3<-c()
for(ff in folders){
  print(ff)
  if(file.exists(paste0(dir,ff,"/plink.ld"))){
    x<-read.table(paste0(dir,ff,"/plink.ld"),header=T, stringsAsFactors = FALSE)
    inv<-unique(as.character(as.matrix(x$SNP_A)))
    x<-x[x$SNP_B!=inv,]
    # x<-x[grep("rs",x$SNP_B),]
    ld<-max(x[order(x$R2,decreasing = T),]$R2)
    
    
    column1<-c(column1,ff)
    column2<-c(column2,ld)
  }else if(length(dir(paste0(dir,ff)))>1){
    coords=read.table("../data/raw/datos_librerias/bplib.coords")
    if(ff=="HsInv0052_region"){
      coords<-c(coords[coords$V1=="REFHsInv0052BP1HG18",3],coords[coords$V1=="REFHsInv0052BP1HG18",4])
    }else{
      coords<-coords[grep(ff,coords$V1),][1,c(3,4)]
    }
    
    
    counter<-1
    for(file in dir(paste0(dir,ff))){
      if(counter==1){
        dfx<-read.table(paste0(dir,ff,"/",file),header=T)
        counter<-counter+1
      }else{
        x<-read.table(paste0(dir,ff,"/",file),header=T)
        dfx<-rbind(dfx,x)
        counter<-counter+1
      }
    }
    
    select_dfx<-unique(dfx[,c(2,3)])
    select_dfx<-as.character(as.matrix(select_dfx[select_dfx$BP_A>(as.numeric(coords[1])-10000) & select_dfx$BP_A<(as.numeric(coords[2])+10000),2]))
    if(ff=="HsInv0052_region"){
      select_dfx<-"BI_GS_CNV_3_162512202_162525833"
    }
    dfx<-dfx[dfx$SNP_A==select_dfx,]
    
    x<-dfx
    inv<-unique(as.character(as.matrix(x$SNP_A)))
    x<-x[x$SNP_B!=inv,]
    x<-x[grep("rs",x$SNP_B),]
    ld<-max(x[order(x$R2,decreasing = T),]$R2)
    
    
    column1<-c(column1,ff)
    column2<-c(column2,ld)
  }
  
  
}
df2<-data.frame(inv=column1,ld_kgp=column2)
df2<-df2[df2$ld_kgp > 0,]
df2<-aggregate(  ld_kgp ~ inv , df2, max)

##### PLOT COMPARISON
not_detected<-setdiff(df$inv,df2$inv)
common<-intersect(df$inv,df2$inv)
rownames(df)<-df$inv
rownames(df2)<-df2$inv
df_c<-df[common,]
df2_c<-df2[common,]
df_c$ld_kgp<-df2_c$ld_kgp

df_c$ld_bs<-1-df_c$ld_bs
df_c$ld_kgp<-1-df_c$ld_kgp

df_plot<-data.frame(inv=c(as.character(as.matrix(df_c$inv)),as.character(as.matrix(df_c$inv))),
                    Approach=c(rep("BS",dim(df_c)[1]),rep("KGP",dim(df_c)[1])),
                    ld=c(df_c$ld_bs,df_c$ld_kgp))
df_plot$inv<-as.character(as.matrix(df_plot$inv))
df_plot$Approach<-as.character(as.matrix(df_plot$Approach))
df_plot$ld<-as.numeric(as.matrix(df_plot$ld))

df_plot<-rbind(df_plot,c("Mean","BS_mean",mean(df_c$ld_bs)))
df_plot<-rbind(df_plot,c("Mean","KGP_mean",mean(df_c$ld_kgp)))
df_plot$Approach<-as.factor(df_plot$Approach)
df_plot$ld<-as.numeric(as.matrix(df_plot$ld))

df_plot_points<-data.frame(inv=df_c$inv,diff=c(df_c$ld_bs-df_c$ld_kgp),Approach="BS")
order<-as.character(as.matrix(df_plot_points[order(df_plot_points$diff),1]))
df_plot$inv<-factor(df_plot$inv,levels=c(order,"Mean"))

ggplot(df_plot, aes(x = as.factor(inv), 
                        y = ld * ((-1)^(Approach == "KGP" | Approach == "KGP_mean")), 
                        fill = Approach)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(name = "Inversion") +
  #geom_point(data=df_plot_points,aes(as.factor(inv),diff),col="black")+
  scale_y_continuous(name = "1-LD", breaks=c(-0.4, -0.2, 0,0.2,0.4),
                     limits = c(-0.5,0.5))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "top")+
  geom_hline(yintercept = 0,linetype=3)+
  scale_fill_manual(values=c("chocolate1","chocolate3", "steelblue2","steelblue4"))
                    

dim(df_plot_points[df_plot_points$diff>0,])
dim(df_plot_points[df_plot_points$diff<0,])
dim(df_plot_points[df_plot_points$diff==0,])


#### Plot not detected
df_plot2<-df[not_detected,]
df_plot2$Approach<-"BS"
colnames(df_plot2)<-c("inv","ld","Approach")
df_plot2<-df_plot2[,c(1,3,2)]
df_plot2$ld<-1-df_plot2$ld

df_plot2$inv<-as.character(as.matrix(df_plot2$inv))
df_plot2$Approach<-as.character(as.matrix(df_plot2$Approach))
df_plot2$ld<-as.numeric(as.matrix(df_plot2$ld))
df_plot2<-rbind(df_plot2,c("Mean","BS_mean",mean(df_plot2$ld)))
df_plot2$Approach<-as.factor(df_plot2$Approach)
df_plot2$ld<-as.numeric(as.matrix(df_plot2$ld))

ggplot(df_plot2, aes(x = as.factor(inv), 
                    y = ld, 
                    fill = Approach)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(name = "Inversion") +
  scale_y_continuous(name = "1-LD", breaks=c(-0.4, -0.2, 0,0.2,0.4),
                     limits = c(-0.5,0.5))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "top")+
  geom_hline(yintercept = 0,linetype=3)+
  scale_fill_manual(values=c("chocolate1","chocolate3", "steelblue2","steelblue4"))+
  annotate("text", x = 7, y = -0.2, label = "Not detected",size=7)


## Boxplot
df_plot3<-df_plot[df_plot$Approach=="BS" | df_plot$Approach=="KGP",]
ggplot(df_plot3, aes(y = Approach, 
                    x = ld,fill=Approach))+ 
  geom_boxplot()+
  scale_x_continuous(name = "1-LD")+
  theme_classic()+
  scale_fill_manual(values=c("chocolate1", "steelblue2"))
  
options(warn=1)
wilcox.test(df_plot3[df_plot3$Approach=="BS",]$ld,df_plot3[df_plot3$Approach=="KGP",]$ld,
            exact = FALSE,alternative = "two.sided")


## Join plot1 and 2
df_plot4<-rbind(df_plot,df_plot2)
df_plot4$Group<-c(rep("Detected",dim(df_plot)[1]),rep("Not detected",dim(df_plot2)[1]))

ggplot(df_plot4, aes(x = (inv), 
                    y = ld * ((-1)^(Approach == "KGP" | Approach == "KGP_mean")), 
                    fill = Approach)) + 
  geom_bar(stat = "identity") +
  scale_x_discrete(name = "Inversion") +
  #geom_point(data=df_plot_points,aes(as.factor(inv),diff),col="black")+
  scale_y_continuous(name = "1-LD", breaks=c(-0.4, -0.2, 0,0.2,0.4),
                     limits = c(-0.5,0.5))+
  theme_classic()+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position = "top")+
  geom_hline(yintercept = 0,linetype=3)+
  scale_fill_manual(values=c("chocolate1","chocolate3", "steelblue2","steelblue4"))+
  facet_grid(. ~ Group, scales = "free_x", space = "free_x")



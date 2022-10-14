setwd("/run/user/1000/gvfs/sftp:host=vedur.hpc.ut.ee,user=jlerga/gpfs/hpchome/jlerga/Imputation/Breakpoints/results/folders_individuals/")
z1<-read.table("z1")
z1<-unique(z1)
inds<-dir(".")
setdiff(inds,unique(as.character(as.matrix(z1$V4)))) # Revisar "V22546" "V49004" : access denied?????
ind<-unique(as.character(as.matrix(z1$V4)))

inversions_test<-c("6","58","59","95","102",
                   "201","379","97","284")
mx<-matrix(,nrow=1,ncol=9)
for(ind_esp in ind){
  z1_esp<-z1[z1$V4==ind_esp,]
  
  #### Ver todos los reads en ambos BPs por conformacion
  z1_esp2<-z1_esp[,c(1,2,4)]
  z1_esp2<-unique(z1_esp2)
  sum_v<-c()
  for(j in 1:dim(z1_esp2)[1]){
    sum_v<-c(sum_v,sum(z1_esp[z1_esp$V1==z1_esp2[j,1] & z1_esp$V2==z1_esp2[j,2] & z1_esp$V4==z1_esp2[j,3],5]))
  }
  z1_esp2<-cbind(z1_esp2,sum_v)
  
  #### Excluir poco fiables (quitar aquellas conformaciones soportadas por menos de X reads)
  z1_excl<-z1_esp2[z1_esp2$sum_v<=3,c(1,2,3)]
  if(dim(z1_excl)[1]==0){
  }else{
    for(k in 1:dim(z1_excl)[1]){
      z1_esp<-z1_esp[!(z1_esp$V1==z1_excl[k,1] & z1_esp$V2==z1_excl[k,2] & z1_esp$V4==z1_excl[k,3]),]
    }
  }
  
  
  #### Alternativa: Excluir según cuantas veces mayor es el numero de reads soportando una conformacion vs contraria
  # support_diff<-20
  # for(k in 1:length(unique(z1_esp2$V1))){
  #   inv_check<-unique(z1_esp2$V1)[k]
  #   
  #   # Medir el support por cada conformacion
  #   inv_support<-z1_esp2[z1_esp2$V1==inv_check & z1_esp2$V2=="INV_REF",4]
  #   std_support<-z1_esp2[z1_esp2$V1==inv_check & z1_esp2$V2=="ST_REF",4]
  #   if(length(std_support)==0 | length(inv_support)==0){
  #   }else{
  #     inv_std<-inv_support/std_support
  #     std_inv<-std_support/inv_support
  #     
  #     # Cogemos cuantas veces mayor es un support respecto del otro
  #     if(inv_std>=std_inv){
  #       support<-inv_std
  #     }else{
  #       support<-std_inv
  #     }
  #     
  #     # Si es mayor que X se quitar la otra conformacion
  #     if(support>support_diff){
  #      if(inv_std>=std_inv){
  #        # quitar reads std como residual
  #        z1_esp<-z1_esp[!(z1_esp$V1==inv_check & z1_esp$V2=="ST_REF"),]
  #      }else{
  #        # quitar reads inv como residual
  #        z1_esp<-z1_esp[!(z1_esp$V1==inv_check & z1_esp$V2=="INV_REF"),]
  #      }
  #     }
  #   }
  # }

  
  
  ##### Hallar genotipos
  genos_v<-c()
  for(inv_esp in inversions_test){
    geno<-""
    genos_ini<-unique(as.character(as.matrix(z1_esp[z1_esp$V1==inv_esp,2])))
    
    if(length(genos_ini)==0){
      geno<-"NA"
    }else if(length(genos_ini)==2){
      geno<-"HET"
    }else{
      if(genos_ini=="ST_REF"){
        geno<-"STD"
      }else{
        geno<-"INV"
      }
    }
    genos_v<-c(genos_v,geno) 
  }
  mx<-rbind(mx,genos_v)
}
mx<-mx[-1,]
rownames(mx)<-ind
colnames(mx)<-c("HsInv0006","HsInv0058","HsInv0059","HsInv0095","HsInv0102",
                   "HsInv0201","HsInv0379","HsInv0097","HsInv0284")
head(mx)
write.table(file="Genotypes_from_BS",mx,col.names = T,row.names = T,quote=F,sep="\t")


##### Comparison
BS_genos<-read.table("Genotypes_from_BS")
tag_genos<-read.table("../Genotypes_from_tagSNPs_preliminar")
common<-intersect(rownames(tag_genos),rownames(BS_genos))
BS_genos<-BS_genos[common,]            
tag_genos<-tag_genos[common,]            

                                                                                  # Quitando    1  2  3  4
all(as.character(as.matrix(BS_genos$HsInv0097))==as.character(as.matrix(tag_genos$HsInv0097)))# T  T  T  T
all(as.character(as.matrix(BS_genos$HsInv0284))==as.character(as.matrix(tag_genos$HsInv0284)))# T  T  T  T
all(as.character(as.matrix(BS_genos$HsInv0379))==as.character(as.matrix(tag_genos$HsInv0379)))# T  T  T  T
all(as.character(as.matrix(BS_genos$HsInv0201))==as.character(as.matrix(tag_genos$HsInv0201)))# F  T  T  T
all(as.character(as.matrix(BS_genos$HsInv0102))==as.character(as.matrix(tag_genos$HsInv0102)))# F  F  T  T
all(as.character(as.matrix(BS_genos$HsInv0095))==as.character(as.matrix(tag_genos$HsInv0095)))# F  F  F  T
all(as.character(as.matrix(BS_genos$HsInv0059))==as.character(as.matrix(tag_genos$HsInv0059)))# F  T  T  T
all(as.character(as.matrix(BS_genos$HsInv0058))==as.character(as.matrix(tag_genos$HsInv0058)))# T  T  T  T
all(as.character(as.matrix(BS_genos$HsInv0006))==as.character(as.matrix(tag_genos$HsInv0006)))# F  F  F  F

# Cuantos fallan?
sum(rep(1,dim(tag_genos)[1])[!as.character(as.matrix(BS_genos$HsInv0201))==as.character(as.matrix(tag_genos$HsInv0201))])
sum(rep(1,dim(tag_genos)[1])[!as.character(as.matrix(BS_genos$HsInv0095))==as.character(as.matrix(tag_genos$HsInv0095))])
sum(rep(1,dim(tag_genos)[1])[!as.character(as.matrix(BS_genos$HsInv0059))==as.character(as.matrix(tag_genos$HsInv0059))])
sum(rep(1,dim(tag_genos)[1])[!as.character(as.matrix(BS_genos$HsInv0058))==as.character(as.matrix(tag_genos$HsInv0058))])
#sum(rep(1,dim(tag_genos)[1])[!as.character(as.matrix(BS_genos$HsInv0006))==as.character(as.matrix(tag_genos$HsInv0006))])
sum(rep(1,dim(tag_genos)[1])[!as.character(as.matrix(BS_genos$HsInv0102))==as.character(as.matrix(tag_genos$HsInv0102))])

sum(rep(1,dim(tag_genos)[1])[!as.character(as.matrix(BS_genos$HsInv0201))==as.character(as.matrix(tag_genos$HsInv0201))])/dim(tag_genos)[1]*100
sum(rep(1,dim(tag_genos)[1])[!as.character(as.matrix(BS_genos$HsInv0095))==as.character(as.matrix(tag_genos$HsInv0095))])/dim(tag_genos)[1]*100
sum(rep(1,dim(tag_genos)[1])[!as.character(as.matrix(BS_genos$HsInv0059))==as.character(as.matrix(tag_genos$HsInv0059))])/dim(tag_genos)[1]*100
sum(rep(1,dim(tag_genos)[1])[!as.character(as.matrix(BS_genos$HsInv0058))==as.character(as.matrix(tag_genos$HsInv0058))])/dim(tag_genos)[1]*100
#sum(rep(1,dim(tag_genos)[1])[!as.character(as.matrix(BS_genos$HsInv0006))==as.character(as.matrix(tag_genos$HsInv0006))])/dim(tag_genos)[1]*100
sum(rep(1,dim(tag_genos)[1])[!as.character(as.matrix(BS_genos$HsInv0102))==as.character(as.matrix(tag_genos$HsInv0102))])/dim(tag_genos)[1]*100

# which are different?
BS_genos[!as.character(as.matrix(BS_genos$HsInv0201))==as.character(as.matrix(tag_genos$HsInv0201)),]
tag_genos[!as.character(as.matrix(BS_genos$HsInv0201))==as.character(as.matrix(tag_genos$HsInv0201)),]

chrom.number<- 1:5
library(tidyverse)
options(scipen=999)
for (i in chrom.number){
#Enter the filtered data and remove rows with missing data
  infilefiltered <- paste("Filter_LD_chromosome_",i,".hap.ld", sep="")
  dataFil <-  read.table(infilefiltered, fill=T, skip=1,skipNul=T)
  dataFil <-  dataFil[,c(2,3,5)]
  dataFil  <- na.omit(dataFil)
#Enter the unfiltered data and remove rows of missing data
  infilenofilter <- paste("NoFilter_LD_chromosome_",i,".hap.ld", sep="")
  dataNoFil <-  read.table(infilenofilter, fill=T, skip=1,skipNul=T)
  dataNoFil <-  dataNoFil[,c(2,3,5)]
  dataNoFil  <- na.omit(dataNoFil)
  
#Transform the data so that it looks like a matrix for filtered and unfiltered data
  LDFil <-pivot_wider(dataFil,names_from = V2, values_from = V5)
  LD2Fil <-LDFil[,-1]
  rownames(LD2Fil) <-c(LDFil$V3)
  matrixLDFil<-data.matrix(LD2Fil)
  
  LDNoFil <-pivot_wider(dataNoFil,names_from = V2, values_from = V5)
  LD2NoFil <-LDNoFil[,-1]
  rownames(LD2NoFil) <-c(LDNoFil$V3)
  matrixLDNoFil<-data.matrix(LD2NoFil)
  
  #Create a Heatmap
  heatmap(matrixLDFil,Colv=NA, Rowv=NA,col=heat.colors(150,rev=T),scale="column",main=paste("Linkage Disequilibrium (Filtered) Chromosome",i))
  dev.copy(png,paste("LD50kb_Heatmap_FilteredData_Chromosome",i,".png",sep=""),width=2000,height=1500,units="px")
  dev.off()
  heatmap(matrixLDNoFil,Colv=NA, Rowv=NA,col=heat.colors(150,rev=T),scale="column",main=paste("Linkage Disequilibrium (Not Filtered) Chromosome",i))
  dev.copy(png,paste("LD50kb_Heatmap_NoFilterData_Chromosome",i,".png",sep=""),width=2000,height=1500,units="px")
  dev.off()
  }


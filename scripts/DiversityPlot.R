options(scipen=999)

#Set window size in bp
windowsize <-100000
chrom.number<-1:5
repeats  <-1:150 

for (i in chrom.number){
  for(x in repeats) { 
 infilesim<- paste("C:/Users/benat/OneDrive/Documents/Biology/Year 4/Dissertation/diversity/diversity_windowof_",windowsize,"_chrm",i,"_Sim",x,".windowed.pi", sep="")
 diversitysim<- read.table(infilesim, fill = T, skip = 1,)
 if (x==1){diversity_table<- diversitysim[,c(2,3,5)]}
  diversity_table[,3+x]<- diversitysim[,5]
  
  if(i==1) {chrm1diversity<-diversity_table}
  if(i==2) {chrm2diversity<-diversity_table}
  if(i==3) {chrm3diversity<-diversity_table}
  if(i==4) {chrm4diversity<-diversity_table}
  if(i==5) {chrm5diversity<-diversity_table}
  }}
diversitydata1<- read.table("C:/Users/benat/OneDrive/Documents/Biology/Year 4/Dissertation/diversity/diversity_windowof_100000_chrm1_Data.windowed.pi", fill = T, skip = 1,)
diversitydata2<- read.table("C:/Users/benat/OneDrive/Documents/Biology/Year 4/Dissertation/diversity/diversity_windowof_100000_chrm2_Data.windowed.pi", fill = T, skip = 1,)
diversitydata3<- read.table("C:/Users/benat/OneDrive/Documents/Biology/Year 4/Dissertation/diversity/diversity_windowof_100000_chrm3_Data.windowed.pi", fill = T, skip = 1,)
diversitydata4<- read.table("C:/Users/benat/OneDrive/Documents/Biology/Year 4/Dissertation/diversity/diversity_windowof_100000_chrm4_Data.windowed.pi", fill = T, skip = 1,)
diversitydata5<- read.table("C:/Users/benat/OneDrive/Documents/Biology/Year 4/Dissertation/diversity/diversity_windowof_100000_chrm5_Data.windowed.pi", fill = T, skip = 1,)

diversitydata1<-diversitydata1[5:301,]
diversitydata2<-diversitydata2[5:191,]
diversitydata3<-diversitydata3[3:231,]
diversitydata4<-diversitydata4[3:185,]
diversitydata5<-diversitydata5[3:259,]

#Obtain the midpoints of each window
diversitydata1$midpoint <- (diversitydata1[[3]]-((windowsize)/2))/1000000
diversitydata2$midpoint <- (diversitydata2[[3]]-((windowsize)/2))/1000000 + max(diversitydata1$midpoint)
diversitydata3$midpoint <- (diversitydata3[[3]]-((windowsize)/2))/1000000 + max(diversitydata2$midpoint)
diversitydata4$midpoint <- (diversitydata4[[3]]-((windowsize)/2))/1000000 + max(diversitydata3$midpoint)
diversitydata5$midpoint <- (diversitydata5[[3]]-((windowsize)/2))/1000000 + max(diversitydata4$midpoint)



diversitydata1$SimMean <-rowMeans(chrm1diversity[,3:ncol(chrm1diversity)])
diversitydata1$StdError <- apply(chrm1diversity[,3:ncol(chrm1diversity)],1,sd,na.rm=T)

diversitydata2$SimMean <-rowMeans(chrm2diversity[,3:ncol(chrm2diversity)])
diversitydata2$StdError <- apply(chrm2diversity[,3:ncol(chrm2diversity)],1,sd,na.rm=T)

diversitydata3$SimMean <-rowMeans(chrm3diversity[,3:ncol(chrm3diversity)])
diversitydata3$StdError <- apply(chrm3diversity[,3:ncol(chrm3diversity)],1,sd,na.rm=T)

diversitydata4$SimMean <-rowMeans(chrm4diversity[,3:ncol(chrm4diversity)])
diversitydata4$StdError <- apply(chrm4diversity[,3:ncol(chrm4diversity)],1,sd,na.rm=T)

diversitydata5$SimMean <-rowMeans(chrm5diversity[,3:ncol(chrm5diversity)])
diversitydata5$StdError <- apply(chrm5diversity[,3:ncol(chrm5diversity)],1,sd,na.rm=T)

diversitydata1$difference <-diversitydata1$V5 - diversitydata1$SimMean

write.table(diversitydata1, "Diversity_windowsof100kbChrm1_Data&Sim.txt", quote=F)
diversitydata2$difference <-diversitydata2$V5 - diversitydata2$SimMean
write.table(diversitydata2, "Diversity_windowsof100kbChrm2_Data&Sim.txt", quote=F)
diversitydata3$difference <-diversitydata3$V5 - diversitydata3$SimMean
write.table(diversitydata3, "Diversity_windowsof100kbChrm3_Data&Sim.txt", quote=F)
diversitydata4$difference <-diversitydata4$V5 - diversitydata4$SimMean
write.table(diversitydata4, "Diversity_windowsof100kbChrm4_Data&Sim.txt", quote=F)
diversitydata5$difference <-diversitydata5$V5 - diversitydata5$SimMean
write.table(diversitydata5, "Diversity_windowsof100kbChrm5_Data&Sim.txt", quote=F)

mean(diversitydata1$V5)
mean(diversitydata1$SimMean)
Totalc1difference <-mean(diversitydata1$V5)-mean(diversitydata1$SimMean)
Totalc2difference <-mean(diversitydata2$V5)-mean(diversitydata2$SimMean)
Totalc3difference <-mean(diversitydata3$V5)-mean(diversitydata3$SimMean)
Totalc4difference <-mean(diversitydata4$V5)-mean(diversitydata4$SimMean)
Totalc5difference <-mean(diversitydata5$V5)-mean(diversitydata5$SimMean)


Diversitydata <- rbind (diversitydata1,diversitydata2,diversitydata3,diversitydata4,diversitydata5)
colnames(Diversitydata) <-c("Chromosome","Bin_Start (bp)", "Bin_End(bp)", "SNP_Number","Real_Diversity","Midpoint","Real_Diversity", "StdError", "Difference")
#Do a statistical analysis of the expected vs observed value in the neutral simulation
t.test(Diversitydata[,7],mu=0.00288176)
mean(Diversitydata[,6]) #2Nemu shoudl be the same as this value. So diversity is as expected.

#Do a statistical test of the diversity in the neutral simulation and the actual data
t.test(Diversitydata[,7],Diversitydata[,5])
#Plot diversity
plot (x=Diversitydata$Midpoint, y=Diversitydata[,5], pch=20, xlab="Mb", ylab=expression(pi),main="Genome-wide diversity in windows of 100kb", type= "l",cex.lab=1.5, cex.axis=1.5, cex.main=2)
lines(x=Diversitydata$Midpoint, y=Diversitydata[,7], pch=20, col="red")


plot (x=diversitydata1$midpoint, y=diversitydata1$V5, pch=20, xlab="Mb", ylab=expression(pi),main="Chromosome 1 diversity in windows of 100kb", type= "l",cex.lab=1.5, cex.axis=1.5, cex.main=2)
lines(x=diversitydata1$midpoint, y=diversitydata1$SimMean, pch=20, col="red")
abline(h= mean(diversitydata1$V5), lty=2)

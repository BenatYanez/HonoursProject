#Do a t-test for each chromosome for the expected and observed ratio of recombination to mutation in the simulations
options(scipen=999)
chrom.number<- 1:5
for (i in chrom.number){
  
  infile <- paste("Chr",i,"_MutationRecombinationRatio.txt", sep="")
  Ratio<- read.table(infile, fill=T, skip=1,skipNul=T)
  LogRatio<- log(Ratio)
  if(i==1){Expected<-0.1228}
  if(i==2){Expected<-0.1645}
  if(i==3){Expected<-0.1457}
  if(i==4){Expected<-0.1715}
  if(i==5){Expected<-0.1369}
  
  LogExpected<-log(Expected)
  #t.test(Ratio,mu=Expected)
  print(paste("T-Test for Cromosome",i))
  print(t.test(Ratio,mu=Expected))
  
}

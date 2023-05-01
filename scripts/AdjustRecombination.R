#Adjust the Salome recombination map to account for the high selfing in Arabidopsis

chrom.number <- 1:5
for (i in chrom.number){
#Open the Salome files and add collumn names
infile <- paste("Salome_Chr",i,"Map_Space.dat", sep="")
Recombination  <- read.table (infile, fill=T, skip=1)
colnames(Recombination) <- c("Chromosome", "pos", "CM")
#Change recombiantion data from CM to recombination probability
Recombination_Prob <- exp(-Recombination$CM/100)*sinh(Recombination$CM/100)
#Adjust recombiantion probability to the expectations of selfing
SelfingRecombination_Prob <- Recombination_Prob*(1-(0.97)/(2-0.97))
#Change the adjusted recombination probability back to CM
CM_Recombination  <- 50*log(1/(1-(2*SelfingRecombination_Prob)))
#Save the data
AdjustedMap<- cbind (Recombination[,1:2], CM_Recombination)
outfile  <- paste("Salome_Chr",i,"Map_Adjusted_Space.dat", sep="")
write.table(AdjustedMap,outfile, quote=F, row.names=F)
}

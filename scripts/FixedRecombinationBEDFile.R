#Create BED files for a fixed recombination size in cM for each chromosome

dat <- read.csv("Salome_AllChrMap.csv")
#Remove scientific notation
options(scipen=999)

#Establish window size in cM
windowsize   <- 2

chrom_number<- 1:5
for (n in chrom_number){
  #Create name of BED file
  outfile <-paste("hapcount_BED_file_chrom",n,"_",windowsize,"cM.bed",sep="")
  cat("track name=\"atsweeps\" description=\"BED file for hapcount chromosome ",n ,"\"\n",
      file = outfile)
#Get the data for just one chromosome
data <- dat[dat$ï..CHROM== n,]
#Reset these values and the BED file everytime you try to run the repeat
cM <- 0
i <- 1
x <- 1


bed   <- data.frame("V1"= c(n),
                    "V2"= 0,
                    "V3"= 488426,
                    "V4"= 0)


data[,4] <- c(0, data[2:nrow(data),3] - data[1:nrow(data)-1,3])
colnames(data)[4] <-"Difference"
repeat {
  
cM <- cM + (data[i,4])

bed[x,3]  <- data[i,2]

i <- i+1
if (is.na(cM) == TRUE ){ break}
#If the Cm are within 25% of windowsize
if ( cM > 1.5 ) {
  bed[x,4] <- cM
  bed <- rbind(bed, c(n,data[i-1,2]+1, data[i,2],0))
  x<- x+1
  cM<- 0
}
}
#Remove entries with NA data, ie, the last one
bed <- na.omit(bed)
bedfile  <- subset(bed, select = 1:3)
write.table(bedfile, file = outfile,
            row.names = F, 
            col.names = F,
            quote = F,
            append = T,)
#Create a file with the intervals, the cM and the mean + standard deviation from it
outtext <-paste("Windows of ",windowsize,"cM size Chromosome",n, " + RecombInfo.txt",sep="")
cat("BED file for hapcount chromosome ",n,"Mean Recombination of windows=",mean(bed[,4]),"Standard Deviation=",sd(bed[,4]),
    file = outtext)
write.table(bed, file = outtext,
            quote= F)

}

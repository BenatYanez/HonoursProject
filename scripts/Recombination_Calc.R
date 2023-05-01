#Calculate recombination rate for each window of windowsize length
dat <- read.csv("Salome_AllChrMap.csv")
#Size of the windows in kb
windowsize  <- 100
chrom.number<-1:5
i<-1
for(i in chrom.number){
  infile <-paste("hapcount_BED_file_chrom",i,"_",windowsize,"kb.bed", sep="")
  outfile <-paste("Rec_AT_Chr",i,"_",windowsize,"kb.dat",sep="")
  #get a dataframe with info for just one chromosome
  data <- dat[dat$ï..CHROM== i,]
  bed <- read.table(infile,sep="",skip=1)
  getRate <- function(map, intI) {
    # get recombination rate in cM/Mb for a given map and interval
    rr1 <- approx( map$pos, map$CM, xout=intI[,1] )
    rr2 <- approx( map$pos, map$CM, xout=intI[,2] )
    rate <- (rr2$y-rr1$y) / ((rr2$x-rr1$x)/1000000) # cM/Mb
    return(rate)
  }
  recout <- cbind(bed[,2:3],getRate(data[,2:3],bed[,2:3]))
  names(recout) <- c("StartPoint","EndPoint","cM/Mb")
  write.table(recout,outfile, quote=F)
  plot(recout[,2]-10000,recout[,3],xlab="Bin Midpoint (Mb)",ylab="cM/Mb"); abline(h=0.5,lty=2)
}


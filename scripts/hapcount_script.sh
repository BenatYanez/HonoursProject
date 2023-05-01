## Calculate H statistics for all rows in the hapcount file using 2 loops.
#Modified from a script by Matthew Fellbaum
# Set working directory
#setwd("~/Documents/Biology/Year4/Dissertation/")

#Set window size in kb
windowsize <-100

chrom.number<-1:5
for(i in chrom.number) {
  #Test if the code is running properly or getting stuck in a specific chromosome. Remove #
  #if(i==1) {print ("Number 1")}
  #if(i==2) {print ("Number 2")}
  #if(i==3) {print ("Number 3")}
  #if(i==4) {print ("Number 4")}
  #if(i==5) {print ("Number 5")}
  hapcount_infile<- paste("hapcount_",windowsize,"kb_chromosome_",i,".hapcount", sep="")
  recomb_infile<- paste("Rec_AT_Chr",i,"_",windowsize,"kb.dat", sep ="")
  
  # Open hapcount table into R, might need to change c(1:14) if there is more collums with hapcount info
  #Depending on window size the number of coloums with haplotype information is different
  hapcount_table <- read.table(hapcount_infile, 
                                 fill = T, skip = 1, col.names = c(1:14),)
#Deal with the extra collumns that fills with NA in hapcounts for chromosomes without 14 collums  
hapcount_table <-replace (hapcount_table, is.na(hapcount_table),"")
  

  # Open recombination table and add recombination column to hapcount table
  rec_table <- read.table(recomb_infile, sep="")
  rec_vector <- rec_table[["cM.Mb"]]
  #Add rec_vector as the extra collumn
  hapcount_table$rec <- rec_vector
  

  #Remove NA data
    hapcount_table <- na.omit(hapcount_table)


  # Filter out rows with SNPs < 100 and cM/Mb < 0.5
  hapcount_table <- hapcount_table[hapcount_table[,4] >= 100,]
  hapcount_table <- hapcount_table[hapcount_table[,ncol(hapcount_table)] >= 0.5,]
  #Remove centromere depending on chromosome
  if (i==1) {
    hapcount_table_1 <- hapcount_table[hapcount_table[,3] <= 13700000,]
    hapcount_table_2 <- hapcount_table[hapcount_table[,2] >= 15900000,]
    hapcount_table <- rbind(x = hapcount_table_1, y = hapcount_table_2)
  }
  if (i==2) {
    hapcount_table_1 <- hapcount_table[hapcount_table[,3] <= 2450000,]
    hapcount_table_2 <- hapcount_table[hapcount_table[,2] >= 5500000,]
    hapcount_table <- rbind(x = hapcount_table_1, y = hapcount_table_2)
  }
  if (i==3) {
    hapcount_table_1 <- hapcount_table[hapcount_table[,3] <= 11300000,]
    hapcount_table_2 <- hapcount_table[hapcount_table[,2] >= 14300000,]
    hapcount_table <- rbind(x = hapcount_table_1, y = hapcount_table_2)
  }
  if (i==4) {
    hapcount_table_1 <- hapcount_table[hapcount_table[,3] <= 1800000,]
    hapcount_table_2 <- hapcount_table[hapcount_table[,2] >= 5150000,]
    hapcount_table <- rbind(x = hapcount_table_1, y = hapcount_table_2)
  }
  if (i==5) {
    hapcount_table_1 <- hapcount_table[hapcount_table[,3] <= 11000000,]
    hapcount_table_2 <- hapcount_table[hapcount_table[,2] >= 13350000,]
    hapcount_table <- rbind(x = hapcount_table_1, y = hapcount_table_2)
  }

    # Define rows and columns
  rows <- nrow(hapcount_table)
  columns <- ncol(hapcount_table)
  row.names(hapcount_table) <- 1:rows
  midpoint <- (hapcount_table[[3]]-((windowsize*1000)/2))/1000000

  # Create hstats vectors
  h1 <- rep(NA, rows)
  h12 <- rep(NA, rows)
  h2 <- rep(NA, rows)
  h21 <- rep(NA, rows)

  # Define loop vectors Makes loop vector 2 have the number 7-14 which correspond to the collums with haplotype info
  loop.vector.1 <- 1:rows
  loop.vector.2 <- 7:(columns-1)


  # Open outer loop for all rows
  for(n in loop.vector.1){
  
    # Define h1 vector for single row to reset in each loop    
    h1.single.row <- rep(NA, length(loop.vector.2))
   
    # Define hapfreq and multiplicity vectors to reset each loop
    hapfreq <- rep(NA, length(loop.vector.2))
    multiplicity <- rep(NA, length(loop.vector.2))
  
    # Open inner loop for 1 row  
    for(f in loop.vector.2){  #For each collumn
    
      # Extract haplotype frequency and multiplicity from hapcount table    
      if(hapcount_table[n , f] == "") {break} #If the If is TRUE then it stops the for function for that row
      r <- hapcount_table[n , f]
      d <- as.numeric(unlist(strsplit(r, ":")))
    
      # Calculate H1 for 1 row    
      j <- ((d[2]/90)^2)*(d[1])  #d[2] is the number of individuals d[1] is the number of haplotypes with that many individuals
      h1.single.row[(f-6)] <- j #stores this calculation in vector
    
      # Define hapfreq and multiplicity vectors for 1 row
      hapfreq[f-6] <- d[2]
      multiplicity[f-6] <- d[1]
     
    } #Done the calculation in every collumn for a single row
  
    # Reverse hapfreq and multiplicity vectors 
    hapfreq <- rev(hapfreq)
    hapfreq <- hapfreq[!is.na(hapfreq)]
    multiplicity <- rev(multiplicity)
    multiplicity <- multiplicity[!is.na(multiplicity)]
  
    # Define 2 most common haplotypes for 1 row 
    hapfreq1 <- 0
    hapfreq2 <- 0
  
    hapfreq1 <- hapfreq[1] 
    if(multiplicity[1] == 1) {hapfreq2 <- hapfreq[2]} 
    #E.g if most common haplotype was 1:16 then to obtain hapfreq2 go to second most common 1:10.
    #If most common haplotype was 2:16 then obtain hapfreq2 agin from the 2:16 since there are two haplotypes with equal number of people
    else {hapfreq2 <- hapfreq[1]}
    if(is.na(hapfreq2)) {hapfreq2 <- 0}
  
    # Calculate H1, H12, H2 and H2/H1 for all rows
    h1[n] <- sum(h1.single.row, na.rm = T)
    h12[n] <- h1[n] + 2*(hapfreq1/90)*(hapfreq2/90) 
    h2[n] <- h1[n]-((hapfreq1/90)^2)
    h21[n] <- h2[n]/h1[n]
  }

  # Bonferroni correction <- 5%/nrows

  bonferroni <- 0.05/rows

  # Calculate quantile of 95% for H values

  h1_quantile <- quantile(h1, probs = 1 - bonferroni)
  h12_quantile <- quantile(h12, probs = 1 - bonferroni)
  h21_quantile <- quantile(h21, probs = 1 - bonferroni)
  
   

  # Identify regions with significant H values.

  h1_table <- cbind(hapcount_table,midpoint, h1, h21)
    #Save tables for each chromosome in the environment
  if(i==1) {chrm1tableh1<-h1_table}
  if(i==2) {chrm2tableh1<-h1_table}
  if(i==3) {chrm3tableh1<-h1_table}
  if(i==4) {chrm4tableh1<-h1_table}
  if(i==5) {chrm5tableh1<-h1_table}

  #Only have significant results in the table
  h1_table <- h1_table[h1_table[, ncol(h1_table)-1] >= h1_quantile,]

  h12_table <- cbind(hapcount_table, h12, h21)
  if(i==1) {chrm1tableh12<-h12_table}
  if(i==2) {chrm2tableh12<-h12_table}
  if(i==3) {chrm3tableh12<-h12_table}
  if(i==4) {chrm4tableh12<-h12_table}
  if(i==5) {chrm5tableh12<-h12_table}
  
  
  h12_table <- h12_table[h12_table[, ncol(h12_table)-1] >= h12_quantile,]

  
  #Save the h1 and h12 tables for each chromosome
  #outfileh1 <-paste("H1_Chrom",i,".txt",sep="")
  #outfileh12 <-paste("H12_Chrom",i,".txt",sep="")
  #write.table(h1_table, outfileh1, quote=F)
  #write.table(h12_table, outfileh12, quote=F)
  
    # Plotting anything above the Bonferroni line can be found in the H1_table and H12_table and then save these as png file
  plottitleh1<-paste("Chromosome ",i," H1 with Bonferroni significance line",sep="")
  plottitleh12<-paste("Chromosome ",i," H12 with Bonferroni significance line",sep="")
  plottitleh21<-paste("Chromosome ",i," H2/H1 with Bonferroni significance line",sep="")
  
  filetittleh1<-paste("Chromosome",i,"_H1_",windowsize,"kb.png",sep="")
  filetittleh12<-paste("Chromosome",i,"_H12_",windowsize,"kb.png",sep="")
  filetittleh21<-paste("Chromosome",i,"_H2-H1_",windowsize,"kb.png",sep="")
  plot(x = midpoint, y = h1, pch = 20, xlab = "Window midpoint (Mb)", main = plottitleh1, ylab = "H1")
  abline(h = h1_quantile, lty = 3)
  dev.copy(png,filetittleh1)
  dev.off()
  
  plot(x = midpoint, y = h12, pch = 20, xlab = "Window midpoint (Mb)", main = plottitleh12, ylab = "H12")
  abline(h = h12_quantile, lty = 3)
  dev.copy(png,filetittleh12)
  dev.off()
  
  plot(x = midpoint, y = h21, pch = 20, xlab = "Window midpoint (Mb)", main = plottitleh21, ylab = "H2/H1")
  abline(h = h21_quantile, lty = 3)
  dev.copy(png,filetittleh21)
  dev.off()
 
}

chrm2tableh1$midpoint <-chrm2tableh1$midpoint + max(chrm1tableh1$midpoint)
chrm3tableh1$midpoint <-chrm3tableh1$midpoint + max(chrm2tableh1$midpoint)
chrm4tableh1$midpoint <-chrm4tableh1$midpoint + max(chrm3tableh1$midpoint)
chrm5tableh1$midpoint <-chrm5tableh1$midpoint + max(chrm4tableh1$midpoint)

genomeh1table <- rbind(chrm1tableh1,chrm2tableh1, chrm3tableh1, chrm4tableh1, chrm5tableh1)
genomeh12table <- rbind(chrm1tableh12,chrm2tableh12, chrm3tableh12, chrm4tableh12, chrm5tableh12)

#Identify significant h1, h12 and h21 values and the top 1% of results
genomerows<- nrow(genomeh1table)
genomebonferroni <- 0.05/genomerows

genomeh1_quantile <- quantile(genomeh1table$h1, probs = 1 - genomebonferroni)
genomeh1_1percent <- genomeh1table[order(genomeh1table$h1, decreasing= TRUE),]
genomeh1_1percent <- head(genomeh1_1percent, genomerows*0.01)

#Save a file with the top 1% of h1 results
write.table(genomeh1_1percent, "Top1%H1values100kb.txt", quote=F)

genomeh12_quantile <- quantile(genomeh12table$h12, probs = 1 - genomebonferroni)
genomeh12_1percent <- genomeh12table[order(genomeh12table$h12, decreasing= TRUE),]
genomeh12_1percent <- head(genomeh12_1percent, genomerows*0.01)
#Save a file with the top 1% of h12 results
write.table(genomeh12_1percent, "Top1%H12values100kb.txt", quote=F)

significant1table <- genomeh1table[genomeh1table[, ncol(genomeh1table)-1] >= genomeh1_quantile,]
significant12table <- genomeh12table[genomeh12table[, ncol(genomeh12table)-1] >= genomeh12_quantile,]
write.table(significant1table, "Significant H1_100kb.txt", quote=F)
write.table(significant12table, "Significant H12_100kb.txt", quote=F)

#Normalize and log transform h1 and h12 statistics 
#Log transform h1, h12 and h2/h1 and normalize it
logh1 <- log(genomeh1table$h1)
normalizedlogh1 <- (logh1- mean(logh1))/sd(logh1)
 #plot(density(normalizedlogh1))
 #qqnorm(normalizedlogh1)
 #qqline(normalizedlogh1)

logh12 <- log(genomeh12table$h12)
normalizedlogh12 <- (logh12- mean(logh12))/sd(logh12)
 #plot(density(normalizedlogh12))
 #qqnorm(normalizedlogh12)
 #qqline(normalizedlogh12)


#Open the file containing the results for the neutral simulations obtained via hstatisticsSimulations script

neutralinfileh1<-paste("NeutralSimulationH1_",windowsize,"kb.txt", sep="")
neutralinfileh12<-paste("NeutralSimulationH12_",windowsize,"kb.txt", sep="")
neutralinfileh21<-paste("NeutralSimulationH21_",windowsize,"kb.txt", sep="")
#Open the results for the neutral simulations
neutralh1 <- read.table(neutralinfileh1, header=T,sep="")
neutralh12<- read.table(neutralinfileh12,header=T,sep="")
neutralh21<- read.table(neutralinfileh21,header=T,sep="")
                                 
#Calculate mean and standard deviations for H-statistic values from neutral simulations
neutralh1$Mean <-rowMeans(neutralh1[,4:(ncol(neutralh1)-1)])
neutralh1$StDev <- apply(neutralh1[,4:ncol(neutralh1)],1,sd,na.rm=T)
neutralh12$Mean <-rowMeans(neutralh12[,4:ncol(neutralh12)])
neutralh12$StDev <- apply(neutralh12[,4:ncol(neutralh12)],1,sd,na.rm=T)
neutralh21$Mean <-rowMeans(neutralh21[,4:ncol(neutralh21)])
neutralh21$StDev <- apply(neutralh21[,4:ncol(neutralh21)],1,sd,na.rm=T)

#Remove rows that have been removed in genome data
#neutralh1  <-neutralh1[-c(134,539,540,541,637,887),]
#neutralh12  <-neutralh12[-c(134,539,540,541,637,887),]
#neutralh21  <-neutralh21[-c(134,539,540,541,637,887),]

 #Information for the creation of custom axis
tickpos  <- vector(mode ="numeric", length(genomerows))
colrs1<- rep("#6E96CE", times = nrow(chrm1tableh1))
colrs2<- rep("#1138AB", times = nrow(chrm2tableh1))
colrs3<- rep("#6E96CE", times = nrow(chrm3tableh1))
colrs4<- rep("#1138AB", times = nrow(chrm4tableh1))
colrs5<- rep("#6E96CE", times = nrow(chrm5tableh1))
colrs <- c(colrs1,colrs2,colrs3,colrs4,colrs5)


tickposmid1 <- max(chrm1tableh1$midpoint)/2
tickposmid2 <- max(chrm1tableh1$midpoint) + ((max(chrm2tableh1$midpoint)-max(chrm1tableh1$midpoint))/2)
tickposmid3 <- max(chrm2tableh1$midpoint) + ((max(chrm3tableh1$midpoint)-max(chrm2tableh1$midpoint))/2)
tickposmid4 <- max(chrm3tableh1$midpoint) + ((max(chrm4tableh1$midpoint)-max(chrm3tableh1$midpoint))/2)
tickposmid5 <- max(chrm4tableh1$midpoint) + ((max(chrm5tableh1$midpoint)-max(chrm4tableh1$midpoint))/2)
tickposmid<- c(tickposmid1,tickposmid2,tickposmid3,tickposmid4,tickposmid5)

#Select the windows containing sweeps from Huber et al., 2015
if (windowsize == 10) {Hubersweeps<- genomeh12table[c(999,1105,1477,3181,4710,5518,5844,7174,8102),]
  Hubersweeps$midpoint<-genomeh1table[c(999,1105,1477,3181,4710,5518,5844,7174,8102),16]}

if (windowsize == 100) {Hubersweeps<- genomeh12table[c(110,124,173,218,372,540,636,674,812,934),]
Hubersweeps$midpoint<-genomeh1table[c(110,124,173,218,372,540,636,674,812,934),16]}

#Plot results
plot(x=genomeh1table$midpoint ,y=genomeh1table$h1,col=colrs, pch = 20, xlab = "Chromosome",xaxt= "n", main = "H1 values for all Chromosomes", ylab = "H1", cex.lab=1.5, cex.axis=1.5, cex.main=2)
  abline(h = genomeh1_quantile, lty = 3)
  lines(x=neutralh1$midpoint, y=neutralh1$Mean, col="red")
  abline(h= min(genomeh1_1percent$h1))
  axis(1,at=tickposmid, labels=labl, cex.axis=1.5)
  dev.copy(png,paste("H1AllChromosomes",windowsize,"kb.png"), width=1000,height=1000, units="px")
  dev.off()
  
plot(x=genomeh1table$midpoint ,y=genomeh12table$h12, col=colrs, pch = 20, xlab = "Chromosome", xaxt= "n", main = "H12 values for all Chromosomes", ylab = "H12", cex.lab=1.5, cex.axis=1.5, cex.main=2)
  abline(h = genomeh12_quantile, lty = 3)
  lines(x=neutralh1$midpoint, y=neutralh12$Mean, col="red")
  abline(h= min(genomeh12_1percent$h12))
  points(x=Hubersweeps$midpoint,y=Hubersweeps$h12,type= "p", pch=25,  col="red", cex=1.5)
  axis(1,at=tickposmid, labels=labl, cex.axis=1.5)
  dev.copy(png,paste("H12AllChromosomes",windowsize,"kb.png"), width=1000,height=800, units="px")
  dev.off()
plot(x=genomeh1table$midpoint ,y=genomeh1table$h21, col=colrs, pch = 20, xlab = "Chromosome",xaxt= "n",  main = "H2/H1 values for all Chromosomes", ylab = "H2/H1", cex.lab=1.5, cex.axis=1.5, cex.main=2)
  axis(1,at=tickposmid, labels=labl, cex.axis=1.5)
  #lines(x=1:nrow(neutralh21), y=neutralh21$Mean, col="red")
  dev.copy(png,paste("H21AllChromosomes",windowsize,"kb.png"), width=1000,height=800, units="px")
  dev.off()
  
  library(scales)

#plot (x=genomeh1table$h1, y=genomeh1table$h21, pch=20, xlab="H1", main="H2/H1 vs H1 values for all Chromosome", ylab="H2/H1", cex.lab=1.5, cex.axis=1.5, cex.main=2)
#points(x=genomeh1_1percent$h1, y=genomeh1_1percent$h21,type= "p", pch=20, col="red")
#dev.copy(png,"H1vsH21.png", width=1000,height=1000, units="px")
#dev.off()

plot (x=genomeh12table$h12, y=genomeh12table$h21, pch=20,cex=2, xlab="H12", main="H2/H1 vs H12 values for all Chromosome", ylab="H2/H1", cex.lab=1.5, cex.axis=1.5, cex.main=2)
points(x=genomeh12_1percent$h12, y=genomeh12_1percent$h21,type= "p",cex=2, pch=20, col="red",)
abline(h=0.35,lty=2)
dev.copy(png,paste("H12vsH21",windowsize,"kb.png"), width=1000,height=800, units="px")
dev.off()

#Model for SNP count and recombination rate vs H1 values

CuratedData<- data.frame(genomeh12table[,c(1,4,15)],normalizedlogh12)
colnames(CuratedData) <- c("Chromosome", "SNPCount", "Recombination", "LogNormalizedH12")
model <- lm(LogNormalizedH12 ~ Recombination + SNPCount + SNPCount*Recombination, data=CuratedData )
plot(model)
summary(model)

plot(x=CuratedData$SNPCount, y=CuratedData$LogNormalizedH12)
plot(x=genomeh12table$X4, y=genomeh12table$h12)

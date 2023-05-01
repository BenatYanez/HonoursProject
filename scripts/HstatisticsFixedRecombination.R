## Calculate H statistics for all rows in the hapcount file using 2 loops for a fixed recombination size.

# Set working directory
#setwd("~/Documents/Biology/Year4/Dissertation/")
#Set window size in cM
windowsize <-2
chrom.number<-1:5
for(i in chrom.number) {
  #Test if the code is running properly or getting stuck in a specific chromosome. Remove #
  if(i==1) {print ("Number 1")}
  if(i==2) {print ("Number 2")}
  if(i==3) {print ("Number 3")}
  if(i==4) {print ("Number 4")}
  if(i==5) {print ("Number 5")}
  hapcount_infile<- paste("hapcount_",windowsize,"cM_chromosome_",i,".hapcount", sep="")
  recomb_infile<- paste("Windows of ",windowsize,"cM size Chromosome",i," + RecombInfo.txt", sep ="")
  
  # Open hapcount table into R, might need to change c(1:14) if there is more collums with hapcount info 
  hapcount_table <- read.table(hapcount_infile, 
                               fill = T, skip = 1, col.names = c(1:11),)
  
  
  # Open recombination table and add recombination column to hapcount table
  rec_table <- read.table(recomb_infile, sep="")
  rec_vector <- rec_table[[4]]

  #Deal with the extra collumn that fills with NA in hapcounts for chromosomes without 14 collums
  hapcount_table <-replace (hapcount_table, is.na(hapcount_table),"")
  #Add rec_vector 
  hapcount_table$rec <- rec_vector
  #Remove NA data- For sample data this removes half the dta since the first half have NA recombination
  hapcount_table <- na.omit(hapcount_table)
  
  
  # Filter out rows with SNPs < 100 and cM/Mb < 0.5
  hapcount_table <- hapcount_table[hapcount_table[,4] >= 100,]
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
  #row.names(hapcount_table) <- 1:rows

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
      if(hapcount_table[n, f] == "") {break} #If the If is TRUE then it stops the for function for that row
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
  
  h1_table <- cbind(hapcount_table, h1, h21)
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
  outfileh1 <-paste("H1_Chrom",i,".txt",sep="")
  outfileh12 <-paste("H12_Chrom",i,".txt",sep="")
  write.table(h1_table, outfileh1, quote=F)
  write.table(h12_table, outfileh12, quote=F)
  
  # Plotting anything above the Bonferroni line can be found in the H1_table and H12_table and then save these as png file
  #plottitleh1<-paste("Chromosome ",i," H1 with Bonferroni significance line",sep="")
  #plottitleh12<-paste("Chromosome ",i," H12 with Bonferroni significance line",sep="")
  #plottitleh21<-paste("Chromosome ",i," H2/H1 with Bonferroni significance line",sep="")
  
  #filetittleh1<-paste("Chromosome",i,"_H1.png",sep="")
  #filetittleh12<-paste("Chromosome",i,"_H12.png",sep="")
  #filetittleh21<-paste("Chromosome",i,"_H2-H1.png",sep="")
  #plot(x = midpoint, y = h1, pch = 20, xlab = "Window midpoint (Mb)", main = plottitleh1, ylab = "H1")
  #abline(h = h1_quantile, lty = 3)
  #dev.copy(png,filetittleh1)
  #dev.off()
  
  #plot(x = midpoint, y = h12, pch = 20, xlab = "Window midpoint (Mb)", main = plottitleh12, ylab = "H12")
  #abline(h = h12_quantile, lty = 3)
  #dev.copy(png,filetittleh12)
  #dev.off()
  
  #plot(x = midpoint, y = h21, pch = 20, xlab = "Window midpoint (Mb)", main = plottitleh21, ylab = "H2/H1")
  #abline(h = h21_quantile, lty = 3)
  #dev.copy(png,filetittleh21)
  #dev.off()
  
  #plot(x = hapcount_table$x15, y = h12, pch = 20, xlab = "Recombination rate (cM/Mb)", main = plottitleRech12, ylab = "H12")
  #dev.copy(png,filetitleRech12)
  #dev.off()
  
  #plot(x = hapcount_table$x15, y = h21, pch = 20, xlab = "Recombination rate (cM/Mb)", main = plottitleRech21, ylab = "H2/H1")
  #dev.copy(png,filetitleRech21)
  #dev.off()
  
}

genomeh1table <- rbind(chrm1tableh1,chrm2tableh1, chrm3tableh1, chrm4tableh1, chrm5tableh1)
genomeh12table <- rbind(chrm1tableh12,chrm2tableh12, chrm3tableh12, chrm4tableh12, chrm5tableh12)

#Identify significant h1, h12 and h21 values and the top 1% of results
genomerows<- nrow(genomeh1table)
genomebonferroni <- 0.05/genomerows

genomeh1_quantile <- quantile(genomeh1table$h1, probs = 1 - genomebonferroni)
genomeh1_1percent <- genomeh1table[order(genomeh1table$h1, decreasing= TRUE),]
genomeh1_1percent <- head(genomeh1_1percent, genomerows*0.01)

write.table(genomeh1_1percent, "Top1%H1values2cM.txt", quote=F)

genomeh12_quantile <- quantile(genomeh12table$h12, probs = 1 - genomebonferroni)
genomeh12_1percent <- genomeh12table[order(genomeh12table$h12, decreasing= TRUE),]
genomeh12_1percent <- head(genomeh12_1percent, genomerows*0.01)
write.table(genomeh12_1percent, "Top1%H12values2cM.txt", quote=F)
significant1table <- genomeh1table[genomeh1table[, ncol(genomeh1table)-1] >= genomeh1_quantile,]
significant12table <- genomeh12table[genomeh12table[, ncol(genomeh12table)-1] >= genomeh12_quantile,]

#Information for the creation of custom axis
tickpos  <- vector(mode ="numeric", length(genomerows))
colrs1<- rep("#6E96CE", times = nrow(chrm1tableh1))
colrs2<- rep("#1138AB", times = nrow(chrm2tableh1))
colrs3<- rep("#6E96CE", times = nrow(chrm3tableh1))
colrs4<- rep("#1138AB", times = nrow(chrm4tableh1))
colrs5<- rep("#6E96CE", times = nrow(chrm5tableh1))
colrs <- c(colrs1,colrs2,colrs3,colrs4,colrs5)

tickpos1 <- nrow(chrm1tableh1)/2
tickpos2 <- length(which(genomeh1table$X1 == 1)) + nrow(chrm2tableh1)/2
tickpos3 <- length(which(genomeh1table$X1 == 1 | genomeh1table$X1 == 2)) + nrow(chrm3tableh1)/2
tickpos4 <- length(which(genomeh1table$X1 == 1 | genomeh1table$X1 == 2 | genomeh1table$X1 == 3)) + nrow(chrm4tableh1)/2
tickpos5 <- length(which(genomeh1table$X1 == 1 | genomeh1table$X1 == 2 | genomeh1table$X1 == 3 | genomeh1table$X1 == 4)) + nrow(chrm5tableh1)/2
tickpos <- c(tickpos1,tickpos2,tickpos3,tickpos4,tickpos5)
labl <- c("1","2","3","4","5")


plot(x=1:nrow(genomeh12table) ,y=genomeh1table$h1,col=colrs, pch = 20,cex=1.5, xlab = "Chromosome",xaxt= "n", main = "H1 values for all Chromosomes", ylab = "H1", cex.lab=1.5, cex.axis=1.5, cex.main=2)
axis(1,at=tickpos, labels=labl, cex.axis=1.5)
dev.copy(png,"H1AllChromosomes.png", width=1000,height=750, units="px")
dev.off()

plot(x=1:nrow(genomeh12table) ,y=genomeh12table$h12, col=colrs, pch = 20, xlab = "Chromosome", xaxt= "n", main = "H12 values for all Chromosomes", ylab = "H12", cex.lab=1.5, cex.axis=1.5, cex.main=2)
abline(h = genomeh12_quantile, lty = 3)
abline(h = min(significantfdr12table$h12), lty = 2)
abline(h= min(genomeh12_1percent$h12))
axis(1,at=tickpos, labels=labl, cex.axis=1.5)
dev.copy(png,"H12AllChromosomes.png", width=1000,height=750, units="px")
dev.off()

plot(x=1:nrow(genomeh1table) ,y=genomeh1table$h21, col=colrs, pch = 20, xlab = "Chromosome",xaxt= "n",  main = "H2/H1 values for all Chromosomes", ylab = "H2/H1", cex.lab=1.5, cex.axis=1.5, cex.main=2)
axis(1,at=tickpos, labels=labl, cex.axis=1.5)
dev.copy(png,"H21AllChromosomes.png", width=1000,height=750, units="px")
dev.off()

plot (x=genomeh1table$h1, y=genomeh1table$h21, pch=20, xlab="H1", main="H2/H1 vs H1 values for all Chromosome", ylab="H2/H1", cex.lab=1.5, cex.axis=1.5, cex.main=2)
dev.copy(png,"H1vsH21.png", width=1000,height=750, units="px")
dev.off()
plot (x=genomeh12table$h12, y=genomeh1table$h21, pch=20, xlab="H12", main="H2/H1 vs H12 values for all Chromosome", ylab="H2/H1", cex.lab=1.5, cex.axis=1.5, cex.main=2)
dev.copy(png,"H12vsH21.png", width=1000,height=750, units="px")
dev.off()



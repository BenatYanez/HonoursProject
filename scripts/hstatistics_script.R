## Calculate H statistics for all rows in the hapcount file using 2 loops.

# Set working directory
#setwd("~/Documents/Biology/Year4/Dissertation/")
chrom.number<-1:5
for(i in chrom.number) {
  #Test if the code is running properly or getting stuck in a specific chromosome. Remove #
  #if(i==1) {print ("Number 1")}
  #if(i==2) {print ("Number 2")}
  #if(i==3) {print ("Number 3")}
  #if(i==4) {print ("Number 4")}
  #if(i==5) {print ("Number 5")}
  hapcount_infile<- paste("hapcount_chromosome_",1,".hapcount", sep="")
  recomb_infile<- paste("Rec_AT_Chr",1,".dat", sep ="")
  
  # Open hapcount table into R, might need to change c(1:14) if there is more collums with hapcount info 
  hapcount_table <- read.table(hapcount_infile, 
                              fill = T, skip = 1, col.names = c(1:14),)
  

  # Open recombination table and add recombination column to hapcount table
  rec_table <- read.table(recomb_infile, sep="")
  rec_vector <- rec_table[["cM.Mb"]]
  #Add rec_vector as the 15th collumn
  hapcount_table$x15 <- rec_vector
  #Deal with the extra collumn that fills with NA in hapcounts for chromosomes without 14 collums
  hapcount_table$X14 <-  replace(hapcount_table$X14,is.na(hapcount_table$X14),"")
  #Remove NA data- For sample data this removes half the dta since the first half have NA recombination
  hapcount_table <- na.omit(hapcount_table)


  # Filter out rows with SNPs < 100 and cM/Mb < 0.5
  hapcount_table <- hapcount_table[hapcount_table[,4] >= 100,]
  hapcount_table <- hapcount_table[hapcount_table[,15] >= 0.5,]
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
  midpoint <- (hapcount_table[[3]]-10000)/1000000

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

  h1_table <- cbind(hapcount_table, h1, h21)
  
  h1_table <- h1_table[h1_table[, 16] >= h1_quantile,]

  h12_table <- cbind(hapcount_table, h12, h21)
  h12_table <- h12_table[h12_table[, 16] >= h12_quantile,]
  
  #Save the h1 and h12 tables for each chromosome
  outfileh1 <-paste("H1_Chrom",i,".txt",sep="")
  outfileh12 <-paste("H12_Chrom",i,".txt",sep="")
  write.table(h1_table, outfileh1, quote=F)
  write.table(h12_table, outfileh12, quote=F)
  
    # Plotting anything above the Bonferroni line can be found in the H1_table and H12_table and then save these as png file
  plottitleh1<-paste("Chromosome ",i," H1 with Bonferroni significance line",sep="")
  plottitleh12<-paste("Chromosome ",i," H12 with Bonferroni significance line",sep="")
  plottitleh21<-paste("Chromosome ",i," H2/H1 with Bonferroni significance line",sep="")
  
  filetittleh1<-paste("Chromosome",i,"_H1.png",sep="")
  filetittleh12<-paste("Chromosome",i,"_H12.png",sep="")
  filetittleh21<-paste("Chromosome",i,"_H2-H1.png",sep="")
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
  
  plottitleRech1<-paste("Chromosome ",i," H1 as a function of recombination rate",sep="")
  plottitleRech12<-paste("Chromosome ",i," H12 as a function of recombination rate",sep="")
  plottitleRech21<-paste("Chromosome ",i," H2/H1 as a function of recombination rate",sep="")
  
  filetitleRech1<-paste("Chromosome",i,"_H1_vs_recombination rate.png",sep="")
  filetitleRech12<-paste("Chromosome",i,"_H12_vs_recombination rate.png",sep="")
  filetitleRech21<-paste("Chromosome",i,"_H2-H1_vs_recombination rate.png",sep="")
 plot(x = hapcount_table$x15, y = h1, pch = 20, xlab = "Recombination rate (cM/Mb)", main = plottitleRech1, ylab = "H1")
  dev.copy(png,filetitleRech1)
  dev.off()
  
  plot(x = hapcount_table$x15, y = h12, pch = 20, xlab = "Recombination rate (cM/Mb)", main = plottitleRech12, ylab = "H12")
  dev.copy(png,filetitleRech12)
  dev.off()
  
  plot(x = hapcount_table$x15, y = h21, pch = 20, xlab = "Recombination rate (cM/Mb)", main = plottitleRech21, ylab = "H2/H1")
  dev.copy(png,filetitleRech21)
  dev.off()
}

#x <- (h1 - mean(h1))/sd(h1)
#plot(density(x))
#qqnorm(x)
#qqline(x)


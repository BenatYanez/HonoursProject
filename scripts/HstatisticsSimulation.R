## Calculate H statistics for all rows in the hapcount file using 2 loops for neutral simulations.

# Set working directory
#setwd("~/Documents/Biology/Year4/Dissertation/")
options(scipen=999)

#Set window size in kb
windowsize <-10
chrom.number<-1:5
#Set the number of replicates obtained from the neutral simulation
repeats  <-1:150

for (i in chrom.number){
  for(x in repeats) { 
  
  #Test if the code is running properly or getting stuck in a specific chromosome. Remove #
  #if(i==1) {print ("Number 1")}
  #if(i==2) {print ("Number 2")}
  #if(i==3) {print ("Number 3")}
  #if(i==4) {print ("Number 4")}
  #if(i==5) {print ("Number 5")}
    hapcount_infile<- paste("C:/Users/benat/OneDrive/Documents/Biology/Year 4/Dissertation/hapcount_sim/hapcount_",windowsize,"kb_chromosome_",i,"_neutralsim_repeat",x,"_Real.hapcount", sep="")
    recomb_infile<- paste("Rec_AT_Chr",i,"_",windowsize,"kbsim.dat", sep ="")
  
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
  #hapcount_table <- hapcount_table[hapcount_table[,4] >= 100,]
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
  
  
  # Identify regions with significant H values.
    if(x==1) {h1_table<- hapcount_table[,1:3] }
    h1_table[,3+x]<-h1
    if(x==1) {h12_table<- hapcount_table[,1:3] }
    h12_table[,3+x]<-h12
    if(x==1) {h21_table<- hapcount_table[,1:3] }
    h21_table[,3+x]<-h21
    #Save tables for each chromosome in the environment
    if(i==1) {chrm1tableh1<-h1_table}
    if(i==2) {chrm2tableh1<-h1_table}
    if(i==3) {chrm3tableh1<-h1_table}
    if(i==4) {chrm4tableh1<-h1_table}
    if(i==5) {chrm5tableh1<-h1_table}
    
    if(i==1) {chrm1tableh12<-h12_table}
    if(i==2) {chrm2tableh12<-h12_table}
    if(i==3) {chrm3tableh12<-h12_table}
    if(i==4) {chrm4tableh12<-h12_table}
    if(i==5) {chrm5tableh12<-h12_table}
        
    if(i==1) {chrm1tableh21<-h21_table}
    if(i==2) {chrm2tableh21<-h21_table}
    if(i==3) {chrm3tableh21<-h21_table}
    if(i==4) {chrm4tableh21<-h21_table}
    if(i==5) {chrm5tableh21<-h21_table}
    
    if(i==1) {chrm1tableh1$midpoint<- midpoint}
    if(i==2) {chrm2tableh1$midpoint<- midpoint}
    if(i==3) {chrm3tableh1$midpoint<- midpoint}
    if(i==4) {chrm4tableh1$midpoint<- midpoint}
    if(i==5) {chrm5tableh1$midpoint<- midpoint}
  
  }
  }

chrm2tableh1$midpoint <-chrm2tableh1$midpoint + max(chrm1tableh1$midpoint)
chrm3tableh1$midpoint <-chrm3tableh1$midpoint + max(chrm2tableh1$midpoint)
chrm4tableh1$midpoint <-chrm4tableh1$midpoint + max(chrm3tableh1$midpoint)
chrm5tableh1$midpoint <-chrm5tableh1$midpoint + max(chrm4tableh1$midpoint)

genomeh1table <- rbind(chrm1tableh1,chrm2tableh1, chrm3tableh1, chrm4tableh1, chrm5tableh1)
genomeh12table <- rbind(chrm1tableh12,chrm2tableh12, chrm3tableh12, chrm4tableh12, chrm5tableh12)
genomeh21table <- rbind(chrm1tableh21,chrm2tableh21, chrm3tableh21, chrm4tableh21, chrm5tableh21)

#Save the neutral simulation results
outfileh1 <-paste("NeutralSimulationH1_",windowsize,"kb.txt",sep="")
write.table(genomeh1table, outfileh1, quote=F)
outfileh12<-paste("NeutralSimulationH12_",windowsize,"kb.txt",sep="")
write.table(genomeh12table, outfileh12, quote=F)
outfileh21<-paste("NeutralSimulationH21_",windowsize,"kb.txt",sep="")
write.table(genomeh21table, outfileh21, quote=F)

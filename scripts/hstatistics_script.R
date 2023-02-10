## Calculate H statistics for all rows in the hapcount file using 2 loops.

# Set working directory
#setwd("~/Documents/Uni/4th_year/Honours_Project/Coding")

# Open hapcount table into R 
hapcount_table <- read.table("haplotype_count_0.98MB.hapcount", 
                             fill = T, skip = 1, col.names = c(1:13))

# Open recombination table and add recombination column to hapcount table
rec_table <- read.table("Rec_AT_098MB.dat")
rec_vector <- rec_table[["cM.Mb"]]
#Add rec_vector as the 14th collumn
hapcount_table$x14 <- rec_vector
#Remove NA data- For sample data this removes half the dta since the first half have NA recombination
hapcount_table <- na.omit(hapcount_table)


# Filter out rows with SNPs < 100 and cM/Mb < 0.5
hapcount_table <- hapcount_table[hapcount_table[,4] >= 100,]
hapcount_table <- hapcount_table[hapcount_table[,14] >= 0.5,]

# Remove the centromere
#hapcount_table_1 <- hapcount_table[hapcount_table[,3] <= 13700000,]
#hapcount_table_2 <- hapcount_table[hapcount_table[,2] >= 15900000,]
#hapcount_table <- rbind(x = hapcount_table_1, y = hapcount_table_2)

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

# Define loop vectors Makes loop vector 2 have the number 7-13 which correspond to the collums with haplotype info
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
h1_table <- h1_table[h1_table[, 15] >= h1_quantile,]

h12_table <- cbind(hapcount_table, h12, h21)
h12_table <- h12_table[h12_table[, 15] >= h12_quantile,]

# Plotting

plot(x = midpoint, y = h1, pch = 20, xlab = "Window midpoint (Mb)", main = "H1 with Bonferroni significance line", ylab = "H1")
abline(h = h1_quantile, lty = 3)

plot(x = midpoint, y = h12, pch = 20, xlab = "Window midpoint (Mb)", main = "H12 with Bonferroni significance line", ylab = "H12")
abline(h = h12_quantile, lty = 3)

plot(x = midpoint, y = h21, pch = 20, xlab = "Window midpoint (Mb)", main = "H2/H1 with Bonferroni significance line", ylab = "H2/H1")
abline(h = h21_quantile, lty = 3)




plot(x = hapcount_table$x14, y = h1, pch = 20, xlab = "Recombination rate (cM/Mb)", main = "H1 as a function of recombination rate", ylab = "H1")
plot(x = hapcount_table$x14, y = h12, pch = 20, xlab = "Recombination rate (cM/Mb)", main = "H12 as a function of recombination rate", ylab = "H12")
plot(x = hapcount_table$x14, y = h21, pch = 20, xlab = "Recombination rate (cM/Mb)", main = "H2/H1 as a function of recombination rate", ylab = "H2/H1")



x <- (h1 - mean(h1))/sd(h1)
plot(density(x))
qqnorm(x)
qqline(x)


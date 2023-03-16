#Create windows of custom length for 5 chromosomes in arabidopsis
#To change window size modify the number in windowsize
#To add more chromosomes copy and paste code from 11-34, change teh chromosome size in lines 13 and 17, 
#change the file in outext1 and change the names of the variables with 1 to the number of the chromosome.

#Remove scientific notation
options(scipen=999)

#Stablish the window size
windowsize  <- 100000

#To use chromosomes of different size chnage 30427613 to the size of the chromosome
chrom1Start <- seq(from = 1, 
                  to = ((floor(30427613/windowsize))*windowsize)-(windowsize-1), 
                  by = windowsize)

chrom1End <- seq(from = windowsize,
                to = (floor(30427613/windowsize))*windowsize,
                by = windowsize)
#Create a column of 1s to put into the table
chrom1 <- rep("1", each = length(chrom1Start))
#Create the bedfile table
bedfile1 <- data.frame(chrom1,chrom1Start,chrom1End)

#Figure out the naming
outext1 <-paste("hapcount_BED_file_chrom1_",windowsize/1000,"kb.bed",sep="")
cat("track name=\"atsweeps\" description=\"BED file for hapcount chromosome 1\"\n",
    file = outext1)

write.table(bedfile1, file = outext1,
            row.names = F, 
            col.names = F,
            quote = F,
            append = T,)

#Create a BED file for Chromosome 2
chrom2Start <- seq(from = 1, 
                   to = ((floor(19698022/windowsize))*windowsize)-(windowsize-1),
                   by = windowsize)

chrom2End <- seq(from = windowsize,
                 to = (floor(19698022/windowsize))*windowsize,
                 by = windowsize)
#Create a column of 2s to put into the table
chrom2 <- rep("2", each = length(chrom2Start))
#Create the bedfile table
bedfile2 <- data.frame(chrom2,chrom2Start,chrom2End)

outext2 <-paste("hapcount_BED_file_chrom2_",windowsize/1000,"kb.bed",sep="")
cat("track name=\"atsweeps\" description=\"BED file for hapcount chromosome 2\"\n",
    file = outext2)

write.table(bedfile2, file = outext2,
            row.names = F, 
            col.names = F,
            quote = F,
            append = T,)
#Create a BED file of chromosome 3
chrom3Start <- seq(from = 1, 
                   to = ((floor(23459780/windowsize))*windowsize)-(windowsize-1),
                   by = windowsize)

chrom3End <- seq(from = windowsize,
                 to = (floor(23459780/windowsize))*windowsize,
                 by = windowsize)
#Create a column of 3s to put into the table
chrom3 <- rep("3", each = length(chrom3Start))
#Create the bedfile table
bedfile3 <- data.frame(chrom3,chrom3Start,chrom3End)

outext3 <-paste("hapcount_BED_file_chrom3_",windowsize/1000,"kb.bed",sep="")
cat("track name=\"atsweeps\" description=\"BED file for hapcount chromosome 3\"\n",
    file = outext3)

write.table(bedfile3, file = outext3,
            row.names = F, 
            col.names = F,
            quote = F,
            append = T,)
#Create a BED file of chromosome 4
chrom4Start <- seq(from = 1, 
                   to = ((floor(18585006/windowsize))*windowsize)-(windowsize-1),
                   by = windowsize)

chrom4End <- seq(from = windowsize,
                 to = (floor(18585006/windowsize))*windowsize,
                 by = windowsize)
#Create a column of 4s to put into the table
chrom4 <- rep("4", each = length(chrom4Start))
#Create the bedfile table
bedfile4 <- data.frame(chrom4,chrom4Start,chrom4End)

outext4 <-paste("hapcount_BED_file_chrom4_",windowsize/1000,"kb.bed",sep="")
cat("track name=\"atsweeps\" description=\"BED file for hapcount chromosome 4\"\n",
    file = outext4)

write.table(bedfile4, file = outext4,
            row.names = F, 
            col.names = F,
            quote = F,
            append = T,)
#Create BED file of chromosome 5
chrom5Start <- seq(from = 1, 
                   to = ((floor(26975439/windowsize))*windowsize)-(windowsize-1),
                   by = windowsize)

chrom5End <- seq(from = windowsize,
                 to = (floor(26975439/windowsize))*windowsize,
                 by = windowsize)
#Create a column of 5s to put into the table
chrom5 <- rep("5", each = length(chrom5Start))
#Create the bedfile table
bedfile5 <- data.frame(chrom5,chrom5Start,chrom5End)

outext5 <-paste("hapcount_BED_file_chrom5_",windowsize/1000,"kb.bed",sep="")
cat("track name=\"atsweeps\" description=\"BED file for hapcount chromosome 5\"\n",
    file = outext5)

write.table(bedfile5, file = outext5,
            row.names = F, 
            col.names = F,
            quote = F,
            append = T,)


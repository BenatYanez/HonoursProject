#Create windows of 20kb, to the length of chromosome 1
chrom1Start <- seq(from = 1, 
                  to = 30400001,
                  by = 20000)

chrom1End <- seq(from = 20000,
                to = 30420000,
                by = 20000)
#Create a column of 1s to put into the table
chrom1 <- rep("1", each = length(chrom1Start))
#Create the bedfile table
bedfile1 <- data.frame(chrom1,chrom1Start,chrom1End)

cat("track name=\"atsweeps\" description=\"BED file for hapcount chromosome 1\"\n",
    file = "hapcount_BED_file_chrom1_20kb.bed")

write.table(bedfile1, file = "hapcount_BED_file_chrom1_20kb.bed",
            row.names = F, 
            col.names = F,
            quote = F,
            append = T,)
options(scipen=999) 

#Create windows of 20kb, to the length of chromosome 2
chrom2Start <- seq(from = 1, 
                   to = 19660001,
                   by = 20000)

chrom2End <- seq(from = 20000,
                 to = 19680000,
                 by = 20000)
#Create a column of 2s to put into the table
chrom2 <- rep("2", each = length(chrom2Start))
#Create the bedfile table
bedfile2 <- data.frame(chrom2,chrom2Start,chrom2End)

cat("track name=\"atsweeps\" description=\"BED file for hapcount chromosome 2\"\n",
    file = "hapcount_BED_file_chrom2_20kb.bed")

write.table(bedfile2, file = "hapcount_BED_file_chrom2_20kb.bed",
            row.names = F, 
            col.names = F,
            quote = F,
            append = T,)
options(scipen=999)
#Create windows of 20kb, to the length of chromosome 3
chrom3Start <- seq(from = 1, 
                   to = 23420001,
                   by = 20000)

chrom3End <- seq(from = 20000,
                 to = 23440000,
                 by = 20000)
#Create a column of 3s to put into the table
chrom3 <- rep("3", each = length(chrom3Start))
#Create the bedfile table
bedfile3 <- data.frame(chrom3,chrom3Start,chrom3End)

cat("track name=\"atsweeps\" description=\"BED file for hapcount chromosome 3\"\n",
    file = "hapcount_BED_file_chrom3_20kb.bed")

write.table(bedfile3, file = "hapcount_BED_file_chrom3_20kb.bed",
            row.names = F, 
            col.names = F,
            quote = F,
            append = T,)
options(scipen=999)
#Create windows of 20kb, to the length of chromosome 4
chrom4Start <- seq(from = 1, 
                   to = 18560001,
                   by = 20000)

chrom4End <- seq(from = 20000,
                 to = 18580000,
                 by = 20000)
#Create a column of 4s to put into the table
chrom4 <- rep("4", each = length(chrom4Start))
#Create the bedfile table
bedfile4 <- data.frame(chrom4,chrom4Start,chrom4End)

cat("track name=\"atsweeps\" description=\"BED file for hapcount chromosome 4\"\n",
    file = "hapcount_BED_file_chrom4_20kb.bed")

write.table(bedfile4, file = "hapcount_BED_file_chrom4_20kb.bed",
            row.names = F, 
            col.names = F,
            quote = F,
            append = T,)
options(scipen=999)
#Create windows of 20kb, to the length of chromosome 5
chrom5Start <- seq(from = 1, 
                   to = 26940001,
                   by = 20000)

chrom5End <- seq(from = 20000,
                 to = 26960000,
                 by = 20000)
#Create a column of 5s to put into the table
chrom5 <- rep("5", each = length(chrom5Start))
#Create the bedfile table
bedfile5 <- data.frame(chrom5,chrom5Start,chrom5End)

cat("track name=\"atsweeps\" description=\"BED file for hapcount chromosome 5\"\n",
    file = "hapcount_BED_file_chrom5_20kb.bed")

write.table(bedfile5, file = "hapcount_BED_file_chrom5_20kb.bed",
            row.names = F, 
            col.names = F,
            quote = F,
            append = T,)
options(scipen=999)
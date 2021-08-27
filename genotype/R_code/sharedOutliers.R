#### Shared outlier SNPs between pairs ####

#Define the path of the outliers per locality
path <- "path-to-file/outlier_SNPs/" 

#Define a vector of pair names (these should be the same as the names in the outlier files)
pairNames <-c("D00H00", "D01H01", "D02H04", "D03H02", "D04H05", "D05H06", "D12H14", "D14H15", "D32H12")

#Create a data frame where rows are all possible outlier SNPs and columns pop pairs, 
#with entries indicating presence/absence of this SNP in the outlier file:

#For each of the population pairs, read the outlier file, and have the SNP in column one, and a 1 in column two, to indicate its presence
for(i in 1:length(pairNames)) {
  SNPs <- data.frame(SNPs = read.table(paste0(path,pairNames[i],"outliers.txt"), header = FALSE),
                     col2 = 1)
  #Add in the column names (SNP and population pair)
  colnames(SNPs) <- c("SNP", pairNames[i])
  if (i==1) {
    allOutlierSNPs <- SNPs
  } else {
    #Merge all files together
    allOutlierSNPs <- merge(allOutlierSNPs, SNPs, all = TRUE)
  }
}
#Replace all NAs with 0
allOutlierSNPs[is.na(allOutlierSNPs)] <- 0
#Sort the table based upon the SNP ID
allOutlierSNPs <- allOutlierSNPs[order(as.vector(allOutlierSNPs$SNP)), ]


#Obtain all pairwise comparisons between pairs, and put this into a data frame with the pairs in the first two columns
pairwiseOutliers <-as.data.frame(t(combn(pairNames, 2)))
#Create four empty columns for storing the number of common SNPs and P-values, common contigs and P-values
pairwiseOutliers <- cbind(pairwiseOutliers,NA,NA)
#Rename the column names
colnames(pairwiseOutliers) <- c("pair1", "pair2", "nCommonSNPs", "pValueSNPs")
#Turn the columns into vectors
pairwiseOutliers$pair1 <- as.vector(pairwiseOutliers$pair1)
pairwiseOutliers$pair2 <- as.vector(pairwiseOutliers$pair2)
#For each of the pairs, calculate the common SNPs, and put the numbers in the pairwiseOutliers table
for(i in 1:nrow(pairwiseOutliers)) {
  #For the SNPs
  pairSNPs <- allOutlierSNPs[,c("SNP", pairwiseOutliers[i,1], pairwiseOutliers[i,2])]
  #Create a new column in the dataframe that multiplies the two columns together. If the result is a 1, this SNP is shared between both pairs
  pairSNPs$common <- pairSNPs[,2] * pairSNPs[,3]
  #Add up the common column i.e. how many common SNPs this pair has. 
  pairwiseOutliers$nCommonSNPs[i] <- sum(pairSNPs$common)
}

#read in the file for n
avSNPsPairs <- read.delim(paste0(path, "avSNPsPairs.txt"), header=FALSE)

#Calculate whether the number of common SNPs is greater than chance. Do this for all pairwise comparisons
#phyper (q, m, n, k...)

for(i in 1:nrow(pairwiseOutliers)) {
  #q = size of overlap minus 1. Found in "commonSNPs"
  q <- pairwiseOutliers$nCommonSNPs[i] - 1
  #m = #outliers in group 1. Counted in each outlier file
  m <- length(readLines(paste0(pairwiseOutliers$pair1[i], "outliers.txt")))
  #n = total # of SNPs minus m (this will be individually calculated for each VCF file, but we could read it from a text document)
  n <- avSNPsPairs$V3[i]-m
  #k = #outliers in group 2. Counted in each outlier file
  k <- length(readLines(paste0(pairwiseOutliers$pair2[i], "outliers.txt")))
  pairwiseOutliers$pValueSNPs[i] <- phyper(q, m, n, k, lower.tail=FALSE)
}
#Note, the number of total SNPs was calculated from the original VCF file which the pairs were extracted from (ESC_rel_50pp_80md_HWE_MAC1_pairs.vcf).
#Doing it this way is a conservative approach. We could also do it with just the number of SNPs considered for each pair
write.table(pairwiseOutliers, file = paste0(path,"sharedSNPsPairwise.txt"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE)


#Per SNP, add up the number of pairs that have it as an outlier
sharedOutliers<-data.frame(SNP=allOutliers$SNPs, Sum=rowSums(allOutliers[2:10]))
#Print this to an output file
write.table(sharedOutliers, file = paste0(path,"sharedSNPsOutliers.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
sharedOutliersCount <- table(sharedOutliers$Sum)




#### Shared outlier genes between pairs ####

#Define the path
path <- "path-to-file/outlier_genes/" 

#Define a vector of pair names 
pairNames <-c("D00H00", "D01H01", "D02H04", "D03H02", "D04H05", "D05H06", "D12H14", "D14H15", "D32H12")

#Create a data frame where rows are all possible outlier genes and columns pop pairs, 
#with entries indicating presence/absence of this gene in the outlier file:

#For each of the population pairs, read the outlier file, and have the gene in column one, and a 1 in column two, to indicate its presence
for(i in 1:length(pairNames)) {
  GENES <- data.frame(SNPs = read.table(paste0(path,pairNames[i],"outlierGenes.txt"), header = FALSE),
                     col2 = 1)
  #Add in the column names (GENE and population pair)
  colnames(GENES) <- c("GENES", pairNames[i])
  if (i==1) {
    allOutlierGENES <- GENES
  } else {
    #Merge all files together
    allOutlierGENES <- merge(allOutlierGENES, GENES, all = TRUE)
  }
}
#Replace all NAs with 0
allOutlierGENES[is.na(allOutlierGENES)] <- 0


#Obtain all pairwise comparisons between pairs, and put this into a data frame with the pairs in the first two columns
pairwiseOutliers <-as.data.frame(t(combn(pairNames, 2)))
#Create four empty columns for storing the number of common genes and P-values, common contigs and P-values
pairwiseOutliers <- cbind(pairwiseOutliers,NA,NA)
#Rename the column names
colnames(pairwiseOutliers) <- c("pair1", "pair2", "nCommonGENES", "pValueGENES")
#Turn the columns into vectors
pairwiseOutliers$pair1 <- as.vector(pairwiseOutliers$pair1)
pairwiseOutliers$pair2 <- as.vector(pairwiseOutliers$pair2)
#For each of the pairs, calculate the common genes, and put the numbers in the pairwiseOutliers table
for(i in 1:nrow(pairwiseOutliers)) {
  #For the SNPs
  pairGENES <- allOutlierGENES[,c("GENES", pairwiseOutliers[i,1], pairwiseOutliers[i,2])]
  #Create a new column in the dataframe that multiplies the two columns together. If the result is a 1, this gene is shared between both pairs
  pairGENES$common <- pairGENES[,2] * pairGENES[,3]
  #Add up the common column i.e. how many common SNPs this pair has. 
  pairwiseOutliers$nCommonGENES[i] <- sum(pairGENES$common)

}

#read in the file for n
avGenesPairs <- read.delim(paste0(path, "avGenesPairs.txt"), header=FALSE)


#Calculate whether the number of common genes is greater than chance. Do this for all pairwise comparisons
#phyper (q, m, n, k...)

for(i in 1:nrow(pairwiseOutliers)) {
  #q = size of overlap minus 1. Found in "commonSNPs"
  q <- pairwiseOutliers$nCommonGENES[i] - 1
  #m = #outliers in group 1. Counted in each outlier file
  m <- length(readLines(paste0(pairwiseOutliers$pair1[i], "outlierGenes.txt")))
  #n = total # of SNPs minus m (this will be individually calculated for each VCF file, but we could read it from a text document)
  n <- avGenesPairs$V3[i]-m
  #k = #outliers in group 2. Counted in each outlier file
  k <- length(readLines(paste0(pairwiseOutliers$pair2[i], "outlierGenes.txt")))
  pairwiseOutliers$pValueGENES[i] <- phyper(q, m, n, k, lower.tail=FALSE)
}
#Note, the number of total genes was calculated from the original VCF file which the pairs were extracted from (ESC_rel_50pp_80md_HWE_MAC1_pairs.vcf).
#Doing it this way is a conservative approach. We could also do it with just the number of genes considered for each pair
write.table(pairwiseOutliers, file = paste0(path,"sharedGenesPairwise.txt"), 
            row.names = FALSE, col.names = TRUE, quote = FALSE)


#Per gene, add up the number of pairs that have it as an outlier
sharedOutliers<-data.frame(GENE=allOutliers$GENES, Sum=rowSums(allOutliers[2:10]))
#Print this to an output file
write.table(sharedOutliers, file = paste0(path,"sharedGenesOutliers.txt"), 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
sharedOutliersCount <- table(sharedOutliers$Sum)


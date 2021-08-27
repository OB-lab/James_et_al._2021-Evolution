library(BaylorEdPsych)

#### Effect size per SNP for ESC pairs ####

#Define the path
path <- "path-to-file/eigens/" 
#Get a list of the input files (all files within the path)
allFiles <- list.files(path = path)
#Create a dataframe ready to store the data
SNPeffectSizes <- data.frame(allFiles)
names <- c("Ecotype", "Pair", "Interaction")
SNPeffectSizes[ , names] <- NA
for(i in 1:length(allFiles)) {
  eigenvectors<-read.delim(paste0(path,allFiles[i]), header=F, sep="")
  head(eigenvectors)
  #In column 1, remove everything except for the population name
  eigenvectors$V1<-gsub("\\-.*","",eigenvectors$V1) 
  #Rename the columns
  colnames(eigenvectors) <- c("Pair", "Individual", "PC1")
  #create the ecotype column
  ecotype<-substr(eigenvectors$Pair, start = 1, stop = 1)
  eigenvectors<-cbind(Ecotype = ecotype, eigenvectors)
  #Replace the D with Dune and H with Headland
  eigenvectors$Ecotype <- gsub('D', 'Dune', eigenvectors$Ecotype)
  eigenvectors$Ecotype <- gsub('H', 'Headland', eigenvectors$Ecotype)
  #Replace the pops with their population pair
  eigenvectors$Pair <- gsub('\\<D00\\>', 'D00H00', eigenvectors$Pair)
  eigenvectors$Pair <- gsub('\\<D01\\>', 'D01H01', eigenvectors$Pair)
  eigenvectors$Pair <- gsub('\\<D02\\>', 'D02H04', eigenvectors$Pair)
  eigenvectors$Pair <- gsub('\\<D03\\>', 'D03H02', eigenvectors$Pair)
  eigenvectors$Pair <- gsub('\\<D04\\>', 'D04H05', eigenvectors$Pair)
  eigenvectors$Pair <- gsub('\\<D05\\>', 'D05H06', eigenvectors$Pair)
  eigenvectors$Pair <- gsub('\\<D12\\>', 'D12H14', eigenvectors$Pair)
  eigenvectors$Pair <- gsub('\\<D14\\>', 'D14H15', eigenvectors$Pair)
  eigenvectors$Pair <- gsub('\\<D32\\>', 'D32H12', eigenvectors$Pair)
  eigenvectors$Pair <- gsub('\\<H00\\>', 'D00H00', eigenvectors$Pair)
  eigenvectors$Pair <- gsub('\\<H01\\>', 'D01H01', eigenvectors$Pair)
  eigenvectors$Pair <- gsub('\\<H04\\>', 'D02H04', eigenvectors$Pair)
  eigenvectors$Pair <- gsub('\\<H02\\>', 'D03H02', eigenvectors$Pair)
  eigenvectors$Pair <- gsub('\\<H05\\>', 'D04H05', eigenvectors$Pair)
  eigenvectors$Pair <- gsub('\\<H06\\>', 'D05H06', eigenvectors$Pair)
  eigenvectors$Pair <- gsub('\\<H14\\>', 'D12H14', eigenvectors$Pair)
  eigenvectors$Pair <- gsub('\\<H15\\>', 'D14H15', eigenvectors$Pair)
  eigenvectors$Pair <- gsub('\\<H12\\>', 'D32H12', eigenvectors$Pair)
  #do the ANOVA
  SNPs <- as.matrix(eigenvectors[,4])
  model <-lm(SNPs~eigenvectors$Ecotype*eigenvectors$Pair)
  #anova(model)
  #extract the effect sizes
  es<-EtaSq(model)
  SNPeffectSizes$Ecotype[i]<-es[4]
  SNPeffectSizes$Pair[i]<-es[5]
  SNPeffectSizes$Interaction[i]<-es[6]
}
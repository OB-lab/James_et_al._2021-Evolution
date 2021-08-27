#### R code for cluster separation score (CSS) ####

#Define the path
path <- "path-to-file/MDS/" 
#Get a list of the input files (all files within the path)
allFiles <- list.files(path = path, pattern = "\\.mds$")
#Create a dataframe to store the data
CSSsnp <- data.frame(allFiles)
names <- c("CSS")
CSSsnp[ , names] <- NA
#Open each file and calculate CSS
for(k in 1:length(allFiles)) {
  MDSall<-read.delim(paste0(path,allFiles[k]), header=T, sep="")
  #Create the ecotype column
  ecotype<-substr(MDSall$IID, start = 1, stop = 1)
  MDSall<-cbind(Ecotype = ecotype, MDSall)
  #Replace the D with Dune and H with Headland
  MDSall$Ecotype <- gsub('D', 'Dune', MDSall$Ecotype)
  MDSall$Ecotype <- gsub('H', 'Headland', MDSall$Ecotype)
  #Subet the data into Dunes and Headlands
  MDSdune <- MDSall[MDSall$Ecotype=='Dune', ]
  MDSheadland <- MDSall[MDSall$Ecotype=='Headland', ]
  
  #Calculate the number of Dune individuals
  m <- nrow(MDSdune)
  #Calculate the number of Headland individuals
  n <-  nrow(MDSheadland)
  
  #Calculate the between group distance
  sumBetween <- 0
  for(i in 1:m) {
    for(j in 1:n) {
      sumBetween <- sumBetween + abs(MDSdune$C1[i] - MDSheadland$C1[j])
    }
  }
  
  #Calculate the within Dune distance
  sumWithinDune <- 0
  for(i in 1:(m-1)) {
    sumWithinDune <- sumWithinDune + abs(MDSdune$C1[i] - MDSdune$C1[i+1])
  }
  
  #Calculate the within Headland distance
  sumWithinHeadland <- 0
  for(j in 1:(n-1)) {
    sumWithinHeadland <- sumWithinHeadland + abs(MDSheadland$C1[j] - MDSheadland$C1[j+1])
  }
  
  #Calculate the CSS
  CSSsnp$CSS[k] <- sumBetween/(m*n) - (m+n)*(sumWithinDune/(m^2 * (m-1)) + sumWithinHeadland/(n^2 * (n-1)))
}

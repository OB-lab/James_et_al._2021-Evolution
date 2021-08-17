#The following R code generates the required files to be used with the CollyerAdamsPCVA function

# 'master.file' is the entire phenotyping file, first column the pair, second column the ecotype,
# third column the population ID, the other columns the phenotype data
master.file <- read.delim ("path-to-file/allTraitsPairs", header = TRUE)

# extract the unique elements of the first column (the pair)
pairs <- unique(master.file$Pair)
  
# create all possible pairwise comparisons
pairwise <-combn(pairs, 2)

# create empty data frame for the results
results <- data.frame(pair1 = pairwise[1,], 
                      pair2 = pairwise[2,],
                      length.v1 = NA,
                      length.v2 = NA,
                      contrast = NA,
                      p.contrast = NA,
                      angle = NA,
                      p.angle = NA)

for(i in 1:ncol(pairwise)) {
  # subset master file to generate proto-y.mat
  y.mat <- master.file[ master.file$Pair %in% pairwise[, i], ]

  # generate x.mat
  x.mat <- matrix(NA, nrow = nrow(y.mat), ncol = 4)  # empty matrix
  x.mat[,1] <- 1  # intercept
  x.mat[,3] <- 2 * (y.mat$Ecotype == "Dune") - 1  # 1 for dunes, -1 for headlands
  x.mat[,2] <- 2 * (y.mat$Pair == pairwise[1,i]) - 1  # 1 for pair #1, -1 for pair #2
  x.mat[,4] <- x.mat[,2] * x.mat[,3]  # interaction
  
  # modify y.mat to fit requirements
  y.mat <- as.matrix(y.mat[, -(1:3)])

  # run CollyerAdamsPCVA function on the two matrices
  results[i, 3:ncol(results)] <- CollyerAdamsPCVA(y.mat, x.mat)
}

write.table(results,file="path-to-file/PCVA_results.txt",row.names=FALSE,col.names=TRUE,quote=FALSE)

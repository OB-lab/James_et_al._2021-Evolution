#### One-way MANOVA ####

# All traits
allTraits <- read.delim ("path-to-file/data_files/allTraits.txt", header = T)
Traits <- as.matrix(allTraits[,4:12])
manovaAll <- manova(Traits ~ allTraits$Ecotype)
manovaAll
summary(manovaAll)

# Plant architecture traits
ArchTraits <- as.matrix(allTraits[,4:8])
manovaArch <- manova(ArchTraits ~ allTraits$Ecotype)
manovaArch
summary(manovaArch)

# Leaf traits
LeafTraits <- as.matrix(allTraits[,9:12])
manovaLeaf <- manova(LeafTraits ~ allTraits$Ecotype)
manovaLeaf
summary(manovaLeaf)


#### Two-way MANOVA ####

allTraitsPairs <- read.delim ("path-to-file/data_files/allTraitsPairs.txt", header = T)
TraitsPairs <- as.matrix(allTraitsPairs[,4:12])
manovaAllPairs <- manova(TraitsPairs ~ allTraitsPairs$Ecotype*allTraitsPairs$Pair)
manovaAllPairs

#Partial effect sizes
library(heplots)
etasq(manovaAllPairs, test="Wilks")


#### K-means clustering ####

# All traits
library(cluster)

allTraits <- read.delim ("path-to-file/data_files/allTraits.txt", header = T)
Traits <- allTraits[4:12]
KmeansPheno <- kmeans(Traits, 2, nstart=25) # k = 2
table(allTraits$Ecotype,KmeansPheno$cluster)

# Plant architecture traits
ArchTraits <- allTraits[4:8]
KmeansPhenoArch <- kmeans(ArchTraits, 2, nstart=25) # k = 2
table(allTraits$Ecotype,KmeansPhenoArch$cluster)

# Leaf traits
LeafTraits <- allTraits[9:12]
KmeansPhenoLeaf <- kmeans(LeafTraits, 2, nstart=25) # k = 2
table(allTraits$Ecotype,KmeansPhenoLeaf$cluster)


#### Linear discriminant analysis  ####

library(MASS)

allTraits <- read.delim ("path-to-file/data_files/allTraits.txt", header = T)
Traits <- as.matrix(allTraits[4:12])
lda(allTraits$Ecotype~Traits) 


#### Vote-counting ####

library(psych)
library(BSDA)

traitMeansPairs <- read.delim ("path-to-file/data_files/traitMeansPairs.txt", header = T)

#Vote-counting for each trait (example is vegetative height):
VH_Dune = traitMeansPairs$VH [traitMeansPairs$Ecotype == "Dune"]
VH_Headland = traitMeansPairs$VH [traitMeansPairs$Ecotype == "Headland"]
test <- SIGN.test(x = VH_Dune, y = VH_Headland, alternative = "two.sided", conf.level = 0.95)
test$p.value


#### Trait-by-trait linear models and effect sizes ####

allTraitsPairs <- read.delim ("path-to-file/data_files/allTraitsPairs.txt", header = T)

#Linear model for each trait (example is vegetative height):
VH_lm <- lm(VH~Ecotype*Pair, data=allTraitsPairs)
anova(VH_lm)

#Effect size for each trait (example is vegetative height):
library(BaylorEdPsych)
EtaSq(VH_lm)


#### Phenotypic Change Vector Analysis ####

#See CollyerAdamsPCVA.R and dataFilesPCVA.R for code


#### Identifying major axes of shared evolutionary change ####

library(rWishart)

nsim <- 100  ## for scampling from Wishart distribution
npop <- 8 ## number of populations
d <- 9 ## number of traits

Xd <- read.csv("path-to-file/AllTraitsTheta.csv", header = FALSE)  ## This the full matrix expanded from the lower off-diagonals the symetric distance matrix (i.e. the pairwise angles between localities)
Xdr< - pi*Xd/180 #convert to radians
C <- cos(Xdr) #convert distance matrix to correlation matrix

eig <- eigen(C)  ### eigen decomposition of C ###
vecs <- eig$vectors  ### eigenvectors of C ###
data <- eig$values  ### eigenvalues of C ###
data_prop <- data/npop ### turning eigenvalues into proportions

#Null expectations, from Wishart distribution
data_null <- matrix(, nrow=nsim, ncol=npop)
r<-rWishart(nsim, d, diag(npop))
for(i in 1:(nsim)){
  C<-cov2cor(r[,,i])
  eig<-eigen(C)
  data_null[i,] <-eig$values
}

#Turning the null expectation eigenvectors into proportions of variance
data_null_prop <- data_null/npop

#Plotting
par(family = "Times New Roman")
plot(data_prop, col = "black", xlab = "Eigenvector", ylab = "Proportion of variance", pch=16, cex.lab=1.3, cex.axis=1.3)
boxplot(data_null_prop,  add = TRUE, boxlty = 1, medlwd = 1, whisklty = 1, medlty = 1, outcex= 0, xaxt = "n", yaxt = "n", boxfill="lightgrey", boxwex=0.3, at = 1:8 + 0)


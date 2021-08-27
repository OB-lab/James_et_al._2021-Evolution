#### Linear models ####

LMs <- read.delim ("path-to-file/phenoEnvGfDivTime.txt", header=TRUE)

#Pheno lengths vs env dist
lm(EnvDist~c, data=LMs)
summary(lm(EnvDist~PhenoLength, data=LMs))
plot(PhenoLength~EnvDist, data=LMs)

#Pheno lengths vs m H2D
lm(mH2Dfw~PhenoLength, data=LMs)
summary(lm(mH2D.fw~PhenoLength, data=LMs))
plot(PhenoLength~mH2D.fw, data=LMs)
plot (lm(mH2Dfw~PhenoLength, data=LMs))

#Pheno lengths vs m D2H
lm(mD2H.fw~PhenoLength, data=LMs)
summary(lm(mD2Hfw~PhenoLength, data=LMs))
plot(PhenoLength~mD2Hfw, data=LMs)

#Pheno lengths vs AVG m
lm(AV.m~PhenoLength, data=LMs)
summary(lm(AVm~PhenoLength, data=LMs))
plot(PhenoLength~AVm, data=LMs)

#Pheno lengths vs divergence time
m <- lm(DivTime~log(PhenoLength), data=LMs)
summary(m)
plot(log(PhenoLength)~DivTime, data=LMs)
plot(m)


#### Mantel tests ####

library(ade4)

# Input data files
SNPs <- read.delim ("path-to-file/sharedSNPsMatrixNoTas.txt", header=FALSE)
genes <- read.delim ("path-to-file/sharedGenesMatrixNoTas.txt", header=FALSE)
pathways <- read.delim ("path-to-file/sharedPathwaysMatrixNoTas.txt", header=FALSE)
angles <- read.delim ("path-to-file/anglesMatrix.txt", header=FALSE)
deltaLengths <- read.delim ("path-to-file/deltaLengthsMatrix.txt", header=FALSE)

SNPs <- as.dist(SNPs)
genes <- as.dist(genes)
pathways <- as.dist(pathways)
angles <- as.dist(angles)
deltaLengths <- as.dist(deltaLengths)

#SNPs vs angles
mantel.rtest(SNPs, angles, nrepet = 999)
# Simulated p-value: 0.986 
plot(angles~SNPs, xlab="Shared outlier SNPs", ylab="Phenotypic angles")
summary(lm(angles~SNPs))
abline(72.01157, -0.16841)

#Genes vs angles
mantel.rtest(genes, angles, nrepet = 999)
# Simulated p-value: 0.886 
plot(angles~genes, xlab="Shared outlier genes", ylab="Phenotypic angles")
summary(lm(angles~genes))
abline(63.40484, -0.15884)

#Pathways vs angles
mantel.rtest(pathways, angles, nrepet = 999)
# Simulated p-value: 0.923
plot(angles~pathways, xlab="Shared enriched pathways", ylab="Phenotypic angles")
summary(lm(angles~pathways))
abline(52.851, -2.703)

#SNPs vs deltaLengths
mantel.rtest(SNPs, deltaLengths, nrepet = 999)
# Simulated p-value: 0.901 
plot(deltaLengths~SNPs, xlab="Shared outlier SNPs", ylab="Phenotypic change in lengths")
summary(lm(deltaLengths~SNPs))
abline(2.282292, -0.005251)

#Genes vs deltaLengths
mantel.rtest(genes, deltaLengths, nrepet = 999)
# Simulated p-value: 0.821
plot(deltaLengths~genes, xlab="Shared outlier genes", ylab="Phenotypic change in lengths")
summary(lm(deltaLengths~genes))
abline(2.072949, -0.005345)

#Pathways vs deltaLengths
mantel.rtest(pathways, deltaLengths, nrepet = 999)
# Simulated p-value: 0.832
plot(deltaLengths~pathways, xlab="Shared outlier pathways", ylab="Phenotypic change in lengths")
summary(lm(deltaLengths~pathways))
abline(1.7761, -0.1028)


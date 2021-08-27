SNPs <-c(1035, 755, 497, 259, 123, 32, 6, 0, 0)
SNPsGenes<-c( 1311, 943, 652, 348, 199, 81, 36, 32, 7)
SNPsFunctions<-c(1040,756,498,262,125,35,6,1,1)
GenesFunctions<-c(281,189,156,92,78,52,30,33,8)
Genes<-c(276,188,155,89,76,49,30,32,7)
Functions<-c(5,1,1,3,2,3,0,1,1)

prop.test(SNPs, SNPsGenes, p = NULL, alternative = "two.sided")
#X2=279.65, df=8, p< 2.2e-16
prop.test(SNPs, SNPsFunctions, p = NULL, alternative = "two.sided")
#X2=361.95, df=8, p< 2.2e-16
prop.test(Genes, GenesFunctions, p = NULL, alternative = "two.sided")
#X2=14.523, df=8, p = 0.06911
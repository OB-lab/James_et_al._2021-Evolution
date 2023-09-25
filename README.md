This repository contains all code and data for the phenotypic and genotypic analysis performed in: **James ME, Wilkinson ME, Bernal DM, Liu H, North HL, Engelstaedter J and Ortiz-Barrientos D. (2021) Phenotypic and genotypic parallel evolution of parapatric ecotypes in *Senecio*. *Evolution*. 75, 3115-3131**

All below code was performed in ```R v3.4.2```, unless otherwise specified. 

# Phenotype

## Leaf morphometrics in Image J

See [ImageJ_leaves.doc](phenotype/protocols/ImageJ_leaves.doc) for the protocol used to extract morphometrics from the leaf samples.

## Datasets

In ~11% of sampled plants, we were unable to measure all six plant architectural traits (such as main stem diameter and main stem angle). In these cases, we took the average of the population to impute the trait value for that individual. (Running the analysis excluding these individuals produced consistent results.)

Each trait was log-transformed and standardised to have a mean of 0 and standard deviation of 1. 

```
dataFrame$trait <- scale(log(dataFrame$trait))
```
Pairwise correlations between all traits were calculated.

```
pairwiseCorr <- cor(traits)
```
, where ```traits``` is a matrix of all phenotypic traits. We removed traits that were highly correlated (>0.8). This data frame of standardised and non-correlated traits was named [allTraits.txt](phenotype/data_files/allTraits.txt)

As some analyses require population pairs, we excluded two Headland allopatric populations (H03 and H07) and two Dune allopatric populations (D09 and D35). This dataset is named [allTraitsPairs.txt](phenotype/data_files/allTraitsPairs.txt).

The vote-counting analyses require the mean trait value per population. For each trait, we calculated the mean trait value for each of the Dune and Headland populations that belong to population pairs. This dataset is named [traitMeansPairs.txt](phenotype/data_files/traitMeansPairs.txt).

See [phenotype.R](phenotype/R_code/phenotype.R) for all phenotypic parallel evolution analyses.

# Genotype

This Genotyping-by-Sequencing dataset comes from previous work (James *et al.* 2021). See https://github.com/OB-lab/James_et_al._2021-MBE for details on DNA extractions, library preparation, sequencing, bioinformatics and data filtering. The vcf file of 9,269 SNPs across all 18 populations was used for all subsequent analyses: [ESC_rel_50pp_80md_HWE_pairs_MAF0.05.vcf](genotype/vcf_file/).

## Linkage disequilibrium

Each population was first extracted into a separate file and we used ```PLINK v1.9``` to calculate the haploblocks per population. For instance, for population D00: 

```
./plink --bfile D03_H02_ESC_rel_50pp_80md_HWE_MAC1_pairs_MAF0.05_D00 --allow-extra-chr --blocks 'no-pheno-req' --out Haploblocks/D00
```

The haploblocks from each population were combined and plotted as a frequency distribution. 

## Parallel nucleotide polymorphisms

To characterize how much genotypic variation of each of the 9,269 sequenced SNPs is explained by the overall differences between ecotypes compared to the individual replicate pairs at each locality, we used ```PLINK``` to conduct a PCA on each SNP. For instance, for SNP tig00000020.18881:

```
./plink --pca 2821 --allow-extra-chr –vcf SNPs/ESC_rel_50pp_80md_HWE_pairs_MAF0.05_tig00000020.18881.vcf --mind –out SNPs/ESC_rel_50pp_80md_HWE_pairs_MAF0.05_tig00000020.18881")
```
The ```.eigenvec``` output files were copied to a separate folder, and we used the R code [SNPeffectSizes.R](genotype/R_code/SNPeffectSizes.R) to calculate the effect size per SNP. These results are here: [SNPeffectSizes.txt](genotype/results/SNP_gene/SNPeffectSizes.txt)

To further explore detailed patterns of parallelism at the level of the nucleotide site, we undertook three approaches. 

#### Approach 1 ####
We detected outliers between the ecotypes (all Dune populations vs all Headland
populations). To do this we undertook three outlier detection methods: the top 1% from the distribution of F<sub>ST</sub> values, the top 1% from the distribution of cluster separation scores (CSS), and those SNPs identified by ```BayeScan```.  

F<sub>ST</sub> per SNP was calculated within ```VCFtools```:
```
vcftools --vcf ESC_rel_50pp_80md_HWE_pairs_MAF0.05.vcf --weir-fst-pop Dune.txt --weir-fst-pop Headland.txt --out ESC_rel_50pp_80md_HWE_pairs_MAF0.05_Dune-Headland
```
, where ```Dune.txt``` and ```Headland.txt``` are lists of the individuals within all the Dune and Headland populations respectively. 

We calculated CSS per SNP by first using ```PLINK``` to perform multidimensional scaling analysis and extracting the first dimension. This was performed separately for each SNP. For instance, for SNP tig00000020.18881:

```
./plink --vcf ESC_rel_50pp_80md_HWE_pairs_MAF0.05_ tig00000020.18881.vcf --cluster --mds-plot 1 --mind --allow-extra-chr –out CSS_SNP/ MDS/ESC_rel_50pp_80md_HWE_pairs_MAF0.05_ tig00000020.18881");
```

We then used a custom R script [CSS.R](genotype/R_code/CSS.R) to calculate the CSS per SNP. 

Outliers were also detected in ```BayeScan``` (the input file was first created using ```PGDspider```, and is found here: [BayescanPairs.txt](genotype/input_files/BayescanPairs.txt)):

```
BayeScan BayeScanPairs.txt -threads 8 pr_odds 10 -o BayeScan_prod100
```
We categorized SNPs as highly differentiated if they contained a posterior probability > 0.91, corresponding to a Bayes Factor of  > 10. For Approach 1, we classified SNPs as outliers if they were detected in at least two of the three methods.

These summarized results from Approach 1 are found here: [DH_SNP_parallelism_summary.xlsx](genotype/results/SNP_gene/SNP_parallelism_summary_DH.xls)

#### Approach 2 ####
We detected outliers separately for the Dune-Headland pairs at each locality (first filtering the VCF file for MAF 0.05 at each locality). As above, we identified outliers for each replicate pair using a combination of F<sub>ST</sub>, CSS, and BayeScan. Instead of selecting the top 1% of F<sub>ST</sub> and CSS values (which is highly dependent on the number of sampled SNPs) we chose stringent cut-offs of F<sub>ST</sub> and CSS > 0.95. Within ```BayeScan```, we considered SNPs with a posterior probability > 0.91 as highly differentiated, which corresponds to a Bayes Factor > 10.

These summarized results from Approach 2 are found here: [SNP_parallelism_summary_DxxHxx.xls](genotype/results/SNP_gene)

#### Approach 3 ####
To detect more subtle signals of outliers between the ecotypes, we asked whether there were concordant allele frequency changes across replicate pairs. Specifically, if a nucleotide site was highly differentiated in at least one pair according to Approach 2, we compared allele frequencies across all pairs for the site. These results are found here: [SNP_concordant_outliers.xls](genotype/results/SNP_gene/SNP_concordant_outliers.xls).

We tested whether the change in allele frequency for each replicate pair was in the same direction across all localities. We used two-sided dependent-samples sign-tests in R to determine the level of statistical significance:   

```
binom.test(9, 9)
# P-value = 0.0039
binom.test(8, 9)
# P-value = 0.039
binom.test(7, 9)
# P-value = 0.1797
```

We therefore treated a SNP as concordant if the change in allele frequency was in the same direction for at least 8 of the 9 localities. 

### *Senecio lautus* transcriptome

To explore whether any of the above candidate outlier SNPs were in genic or non-genic regions, we created the *S. lautus* transcriptome. The raw RNAseq files are available upon request. We first used ```SEECER v0.1.3``` with default parameters to correct errors in RNAseq data.

```
bash ./bin/run_seecer.sh reads1.fastq reads2.fastq 
```

We then carried out de novo assembly using ```Trinity v2.0.2``` :

```
perl /path/to/Trinity.pl --seqType fq --left R1.fastq --right R2.fastq --JM 2G --min_contig_length 100 --CPU 4 --bflyHeapSpaceMax 10G
```
, where ```JM``` determines the amount of memory used for k-mer counting by jellyfish, and ```bflyHeapSpaceMax``` defines the amount of memory that is allocated to a single Butterfly (Trinity compoment) process. 

We chose the representative transcripts from each locus by retaining the transcript with the highest read coverage for each subcomponent. To further remove the redundant transcripts, which may come from alternative splicing or close paralogs, we clustered the assembly using ```CD-HIT-EST v4.6.1``` 

```
cd-hit-est -i in.fasta -o out.fasta -c 96 -n 8 -r 1 
```
, where ```c``` is the sequence identity threshold, ```n``` is the word length, and ```r``` is a value of either  1 or 0, where 1 represents that the chosen alignment will be done on both + and - strands.

We then chose the longest transcript as a representative from each cluster.

This assembled transcriptome was mapped to the *S. lautus* PacBio reference genome using ```minimap2```:

```
minimap2 Senecio_PacBio_v1.fasta Senecio_transcriptome_v1.fasta > mapped.paf
```

A series of scripts and file manipulations (details available upon request) were undertaken to determine whether SNPs resided in genic or non-genic regions. More specifically, we considered each transcript a separate gene, which included all isoforms. As the transcriptome excludes introns, we still considered SNPs mapped to the reference genome that fall between two segments of the same transcript as a genic SNP. All other SNPs were considered non-genic, which are expected to include variants in regulatory and repetitive regions as well as in genic regions with unknown homologous genes in other plants. We excluded SNPs that had > 1 gene mapping to it.

We calculated the number of shared outlier SNPs between all pairwise comparisons across localities using the “Shared outlier SNPs between pairs” section of the R script: [sharedOutliers.R](genotype/R_code/sharedOutliers.R). This script also calculates whether the number of shared outliner SNPs is greater than chance. The input files (i.e., the outlier SNPs per locality, “DxxHxxoutliers.txt”) are found here: [outlier_SNPs](genotype/input_files/outlier_SNPs). Also in this folder is the [avSNPsPairs.txt](genotype/input_files/outlier_SNPs/avSNPsPairs.txt) file, which was used for the total number of SNPs per locality (see R code for details). These results are found here: [sharedSNPsPairwise.txt](genotype/results/SNP_gene/sharedSNPsPairwise.txt)

We also calculated the number of total pairs that have a given SNP as an outlier (see R code for details). These results are found here: [sharedSNPsOutliers.txt](genotype/results/SNP_gene/sharedSNPsOutliers.txt)

## Parallel genic polymorphisms

As with the SNPs above, we characterized how much genotypic variation of each of the genes is explained by the overall differences between ecotypes compared to the individual replicate pairs at each locality, we used ```PLINK``` to conduct a PCA on each gene. For instance, for gene comp66_c0.vcf:

```
./plink --pca 2821 --allow-extra-chr –vcf genes/comp66_c0.vcf --mind –out genes/ESC_rel_50pp_80md_HWE_pairs_MAF0.05_tig00000020.18881")
```
, where ```comp66_c0.vcf``` is a VCF file containing all the sequenced SNPs in that gene. We retained the loadings of the first eigenvector for each gene.

The ```.eigenvec``` output files were copied to a separate folder, and we used the same R code as specified above [SNPeffectSizes.R](genotype/R_code/SNPeffectSizes.R) to calculate the effect size per gene. These results are here: [geneEffectSizes.txt](genotype/results/SNP_gene/geneEffectSizes.txt)

For each outlier gene of interest, we obtained it’s RefSeq code by using ```BLASTx``` with the *S. lautus* transcript in which the outlier SNP fell within. 

```
blastx -query genes_transcripts.fasta -db refseq_protein -entrez_query "Arabidopsis thaliana [organism]" -remote -evalue 1e-6 -outfmt 7 -out genesAnnoAll.txt
```
, where ```genes_transcripts.fasta``` are the *S. lautus* transcripts of interest. We then extracted the top hit for each transcript. These results are here: [genesTopHitsAll.txt](genotype/results/SNP_gene/genesTopHitsAll.txt)

We used ```DAVID``` to obtain the predicted functional annotation for the outlier genes of interest. These are found in the “DxxHxx-ClusteringSummary.xls” files, which are explained below. 

We calculated the number of shared outlier genes between all pairwise comparisons across localities using the “Shared outlier genes between pairs” section of the R script: [sharedOutliers.R](genotype/R_code/sharedOutliers.R). This script also calculates whether the number of shared outliner genes is greater than chance. 

The input files (i.e., the outlier genes per locality, “DxxHxxoutlierGenes.txt”) are found here: [outlier_genes](genotype/input_files/outlier_genes). Also in this folder is the following file [avGenesPairs.txt](genotype/input_files/outlier_genes/avGenesPairs.txt), which is used for the total number of genes per locality (see R code for details). These results are found here: [sharedGenesPairwise.txt](genotype/results/SNP_gene/sharedGenesPairwise.txt)

We also calculated the number of total pairs that have a given gene as an outlier (see R code for details). These results are found here: [sharedGenesOutliers.txt](genotype/results/SNP_gene/sharedGenesOutliers.txt)

## Enriched biological functions

We undertook gene-enrichment analysis for the outlier genes for each replicate pair using functional annotation clustering in the web-based program ```DAVID```, using the *Arabidopsis* orthologues for our outlier genes. We used the *Arabidopsis thaliana* genome as the genetic background. These input files per pair are found here: [genotype/input_files/functional_pathway](genotype/input_files/functional_pathway/). 

The results per pair are found here: [genotype/results/functional_pathway](genotype/results/functional_pathway/). The first sheet of each excel document (“DxxHxx-ClusteringSummary.xls”) is the direct output from ```DAVID```. In the P-value column, cells are highlighted green if the P-value < 0.05. The second sheet “TopFromEachCluster” contains the category with the smallest P-value for each cluster (i.e., the rows in bold in the first sheet of the spreadsheet). An additional column “Summarised term” was added to summarise the terms from each cluster when there were multiple P-values < 0.05 within a cluster. 

[allPairs-ClusteringSummary.xls](genotype/results/functional_pathway/allPairs-ClusteringSummary.xls) shows how these summarised terms were grouped into their final functional pathway categories.

## Distributions of shared SNPs, genes, pathways

We compared the distributions of the proportions of shared outlier nucleotide sites, outlier genes and enriched biological functions across pairs using a two-sided X<sub>2</sub>-test with continuity correction in R: [chisq.R](genotype/R_code/chisq.R).
See: [PropSharedPairs.txt](genotype/results/SNP_gene_pathway/PropSharedPairs.txt) for the summary file of the counts and proportions of shared SNPs, genes and functions across pairs. 

## Demographic effects on phenotypic parallelism

We tested whether the variation in phenotypic parallelism within the system could be explained by demographic factors. See the section “Linear models” of the R code: [phenotypeGenotype.R](phenotype_genotype/R_code/phenotypeGenotype.R) for the linear models. The input file is found here: [phenoEnvGfDivTime.txt](phenotype_genotype/input_files/phenoEnvGfDivTime.txt) (note some of these values are found in Table S4 of the manuscript). 
We asked whether pairs that were more phenotypically similar shared more outlier nucleotide sites, genes, and biological functions using Mantel tests. The input files are as follows:
The matrix of phenotypic angles between localities: [anglesMatrix.txt](phenotype_genotype/input_files/anglesMatrix.txt)

The matrix of phenotypic change in lengths between localities: [deltaLengthsMatrix.txt](phenotype_genotype/input_files/deltaLengthsMatrix.txt)

The matrix of shared outlier SNPs between localities: [sharedSNPsMatrixNoTas.txt](phenotype_genotype/input_files/sharedSNPsMatrixNoTas.txt)

The matrix of shared outlier genes between localities: [sharedGenesMatrixNoTas.txt](phenotype_genotype/input_files/sharedGenesMatrixNoTas.txt)

The matrix of shared enriched biological functions between localities: [sharedPathwaysMatrixNoTas.txt](phenotype_genotype/input_files/sharedPathwaysMatrixNoTas.txt)

In these files, the populations are ordered numerically (i.e., D00H00, D01H01, D02H04, D03H02, D04H05, D05H06, D12H14, D14H15, D32H12). See the section “Mantel tests” of the R code: [phenotypeGenotype.R](phenotype_genotype/R_code/phenotypeGenotype.R) for the Mantel tests. 










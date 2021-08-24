This repository contains all code and data for the phenotypic and genotypic analysis performed in: **James ME, Wilkinson ME, Bernal DM, Liu H, North HL, Engelstaedter J and Ortiz-Barrientos D. (2021) Phenotypic and genotypic parallel evolution of parapatric ecotypes in *Senecio*. [journal] [volume.page] [doi]**

All below code was performed in ```R v3.4.2```, unless otherwise specified. 

# Phenotype

## Leaf morphometrics in Image J

See [ImageJ_leaves.doc](phenotype/protocols/ImageJ_leaves.doc) for the protocol used to extract morphometrics from the leaf samples.

## Datasets

In ~11% of sampled plants, we were unable to measure all six plant architectural traits (such as main stem diameter and main stem angle). In these cases, we took the average of the population to impute the trait value for that individual. (Running the analysis excluding these individuals produced consistent results.)

Each trait was log-transformed and standardised to have a mean of 0 and standard deviation of 1. 

```
dataFrame$trait <- scale(log(dataFrame$trait)
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

This Genotyping-by-Sequencing dataset comes from previous work (James *et al.* 2021). See https://github.com/OB-lab/James_et_al._2021-MBE for details on DNA extractions, library preparation, sequencing, bioinformatics and data filtering. The vcf file of 9,269 SNPs across all 18 populations was used for all subsequent analyses: [ESC_rel_50pp_80md_HWE_pairs_MAF0.05.vcf](genotype/vcf_file/ESC_rel_50pp_80md_HWE_pairs_MAF0.05.vcf).

## Linkage disequilibrium

Each population was first extracted into a separate file, and we used ```PLINK v1.9``` to calculate the haploblocks per population. For instance, for population D00: 

```
./plink --bfile D03_H02_ESC_rel_50pp_80md_HWE_MAC1_pairs_MAF0.05_D00 --allow-extra-chr --blocks 'no-pheno-req' --out Haploblocks/D00
```

The haploblocks from each population were combined and plotted as a frequency distribution. 

## Parallel nucleotide polymorphisms

To characterize how much genotypic variation of each of the 9,269 sequenced SNPs is explained by the overall differences between ecotypes compared to the individual replicate pairs at each locality, we used ```PLINK``` to conduct a PCA on each SNP. For instance, for SNP tig00000020.18881:

```
./plink --pca 2821 --allow-extra-chr –vcf SNPs/ESC_rel_50pp_80md_HWE_pairs_MAF0.05_tig00000020.18881.vcf --mind –out SNPs/ESC_rel_50pp_80md_HWE_pairs_MAF0.05_tig00000020.18881")
```
The ```.eigenvec``` output files were copied to a separate folder, and we used the R code [SNPeffectSizes.R](genotype/R_code/SNPeffectSizes.R) to calculate the effect size per SNP.  These results are here: [SNPeffectSizes.txt](genotype/output_files/SNPeffectSizes.txt)

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

We then used a custom R script [CSS.R](R_code/CSS.R) to calculate the CSS per SNP. 

Outliers were also detected in ```BayeScan``` (the input file was first created using ```PGDspider```):

```
BayeScan ESC_rel_50pp_80md_HWE_pairs_MAF0.05.BayeScan.txt -threads 8 pr_odds 10 -o BayeScan_prod100
```
For ```Approach 1```, we classified SNPs as outliers if they were detected in at least two of the three methods.

These summarized results from ```Approach 1``` are found here: [DH_SNP_parallelism_summary.xlsx](genotype/results/DH_SNP_parallelism_summary.xls)

#### Approach 2 ####
We detected outliers separately for the Dune-Headland pairs at each locality. As above, we identified outliers for each replicate pair using a combination of F<sub>ST</sub>, CSS, and BayeScan. Instead of selecting the top 1% of F<sub>ST</sub> and CSS values (which is highly dependent on the number of sampled SNPs) we chose stringent cut-offs of F<sub>ST</sub> and CSS > 0.95. Within ```BayeScan```, we considered SNPs with a posterior probability > 0.91 as highly differentiated, which corresponds to a Bayes Factor > 10.

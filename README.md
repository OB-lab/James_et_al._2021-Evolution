This repository contains all code and data for the phenotypic and genotypic analysis performed in: **James ME, Wilkinson ME, Bernal DM, Liu H, North HL, Engelstaedter J and Ortiz-Barrientos D. (2021) Phenotypic and genotypic parallel evolution of parapatric ecotypes in *Senecio*. [journal] [volume.page] [doi]**

All below code was performed in ```R v3.4.2```, unless otherwise specified. 

# Phenotype

## Leaf morphometrics in Image J

See [ImageJ_leaves.doc](protocols/ImageJ_leaves.doc) for the protocol used to extract morphometrics from the leaf samples.

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
, where ```traits``` is a matrix of all phenotypic traits. We removed traits that were highly correlated (>0.8). This data frame of standardised and non-correlated traits was named [allTraits.txt](data_files/allTraits.txt)

As some analyses require population pairs, we excluded two Headland allopatric populations (H03 and H07) and two Dune allopatric populations (D09 and D35). This dataset is named [allTraitsPairs.txt](data_files/allTraitsPairs.txt).

The vote-counting analyses require the mean trait value per population. For each trait, we calculated the mean trait value for each of the Dune and Headland populations that belong to population pairs. This dataset is named [traitMeansPairs.txt](data_files/traitMeansPairs.txt).

See [phenotype.R](R_code/phenotype.R) for all phenotypic parallel evolution analyses.

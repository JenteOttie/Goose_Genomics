# Workflow for the Bean Goose Project

This project investigates the evolutionary history of Taiga Bean Goose (*Anser fabalis*) and Tundra Bean Goose (*Anser serrirostris*). The dataset is comprised of 9 individuals for each species.

## Mapping to reference genome and SNP calling

This pipeline is described in the [Data_Processing](https://github.com/JenteOttie/Goose_Genomics/blob/master/Data_Processing.md) section.

## Population Structure: PCA and Admixture

To assess the population structure of the Bean Geese, I will make PCAs and Admixture plots using different settings to filter the SNPs generated in the previous step. The filtering will be done using Plink version 1.07 and VCFtools version 0.1.15.

The worflow for generating the PCA can be found [here](https://github.com/JenteOttie/Goose_Genomics/blob/master/BeanGoose/PCA.md)

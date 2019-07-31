# Workflow for the Bean Goose Project

This project investigates the evolutionary history of Taiga Bean Goose (*Anser fabalis*) and Tundra Bean Goose (*Anser serrirostris*). The dataset is comprised of 9 individuals for each species.

&nbsp;

## Mapping to reference genome and SNP calling

This pipeline is described in the [Data_Processing](https://github.com/JenteOttie/Goose_Genomics/blob/master/Data_Processing.md) section.

&nbsp;

## Population Structure: PCA and Admixture

To assess the population structure of the Bean Geese, I will make PCAs and Admixture plots using different settings to filter the SNPs generated in the previous step. The filtering will be done using Plink version 1.07 and VCFtools version 0.1.15.

The workflow for generating the PCA-plots can be found [here](https://github.com/JenteOttie/Goose_Genomics/blob/master/BeanGoose/PCA.md). 
Changing the filtering parameters does not markedly influence the patterns in the PCA. With a more strict filtering (and thus less SNPs), the two species remain distinct. However, there is less structure within the species after strict filtering.

The workflow for generating the Admixture-plots can be found here.

&nbsp;

## Window Analyses

Summary statistics calculated in windows of different sizes (10kb, 20kb, 200kb).

&nbsp;

## The Genomic landscape

Visualising the genomic landscape based on 200kb windows.

&nbsp;

## Demographic analyses

Reconstructing the evolutionary history using dadi.

&nbsp;

## Finding introgression regions using Machine Learning

In the previous step - demographic analysis - I showed that Taiga and Tundra Bean Goose have exchanged DNA during their evolutionary history. To explore which genomic regions have been exchanged, I will apply a machine learning approach: FILET (Finding Introgressed Loci via Extra-Trees). The original paper can be found [here](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007341).

The workflow for the machine learning analysis can be found [here](https://github.com/JenteOttie/Goose_Genomics/blob/master/BeanGoose/Machine%20Learning.md).

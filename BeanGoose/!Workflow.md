# Workflow for the Bean Goose Project

This project investigates the evolutionary history of Taiga Bean Goose (*Anser fabalis*) and Tundra Bean Goose (*Anser serrirostris*). The dataset is comprised of 9 individuals for each species.

## Mapping to reference genome and SNP calling

This pipeline is described in the [Data_Processing](https://github.com/JenteOttie/Goose_Genomics/blob/master/Data_Processing.md) section.

## Population Structure: PCA and Admixture

To assess the population structure of the Bean Geese, I will make PCAs and Admixture plots using different settings to filter the SNPs generated in the previous step. The filtering will be done using Plink version 1.07 and VCFtools version 0.1.15. The steps below show one example with settings for Hardy-Weinberg, MAF and LD.

**Convert VCF-file to Plink-format**

$ vcftools --vcf BeanGoose.SNPs.HardFilt.vcf --plink --out BeanGoose

**Make binary ped-files (bed)**

$ plink –file BeanGoose –make-bed –out BeanGoose

**Filter file on Hardy-Weinberg (0.01) and minor allele frequencies (0.05)**

$ plink --bfile BeanGoose --hwe 0.01 --maf 0.05 --make-bed --out BeanGooseFiltered

**Prune samples based on linkage disequilibrium (sliding window of 50 SNPs, shifting 10 at a time with a R2 threshold of 0.05)**

Because Plink cannot work with too many scaffolds, I will first change all scaffolds to 1.

$ awk '$1="1"' $file\HWE_MAF.bim >temp.bim # Change scaffold column to 1

$ rm $file_HWE_MAF\.bim # Remove old bim-file

$ mv temp.bim $file_HWE_MAF\.bim # Change name of temporary bim-file

Now everything is read to calculate LD for all the SNPs and prune them

$ plink --bfile BeanGooseFiltered --indep-pairwise 50 10 0.05 --make-bed --out BeanGooseFiltered

$ plink --bfile BeanGooseFiltered --extract $BeanGooseFiltered.in --make-bed --out BeanGooseFiltered

**Make PCA**

$ plink --bfile BeanGooseFiltered --pca 4 --out BeanGooseFiltered

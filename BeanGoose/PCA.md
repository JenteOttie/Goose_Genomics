# Generating PCA-plots with different filtering options
Based on the script PCA_generator.sh
## Filtering SNPs
Variants will be filtered based on different criteria:
-Deviations from Hardy-Weinberg equilibrium (HWE)
-Minor allele frequency (MAF)
-Linkage disequilibrium (LD)
The thresholds for these filter can be adjusted.

&nbsp;

**Set parameters for filtering**
```
HWE_parameter=0.01
MAF_parameter=0.05
LD_parameter=0.5
```

**Load input-file**
The input-file is a hard-filtered VCF-file. 

The filename should end on '..SNPs.HardFilt.vcf'.
```
file=$(basename $1 .SNPs.HardFilt.vcf)
echo 'Loading file'
```

**Convert VCF-file to Plink-format**

The VCF-file is converted into Plink-format, only keeping biallelic SNPs.
```
vcftools --vcf $file\.SNPs.HardFilt.vcf --plink --out $file
echo 'Converted file to Plink-format'
```

**Make BED-file**
```
plink --file $file --make-bed --out $file
echo 'Made BED-file'
```

**Filter file on Hardy-Weinberg and minor allele frequencies (set parameters above)**
```
plink --bfile $file --hwe $HWE_parameter --maf $MAF_parameter --make-bed --out $file\_HWE_MAF
echo 'Filtered om Hardy-Weinberg and MAF'
```

**Prune samples based on linkage disequilibrium (LD)**

First, set all chromosomes to 1 in bim-file - Plink cannot work with too many scaffolds
```
awk '$1="1"' $file\_HWE_MAF.bim >temp.bim # Change scaffold column to 1
rm $file\_HWE_MAF.bim # Remove old bim-file
mv temp.bim $file\_HWE_MAF.bim # Change name of temporary bim-file
```

**Calculate LD for samples in sliding windows of 50 SNPs using steps of 10 SNPs. LD-parameter can be changed above**
```
plink --bfile $file\_HWE_MAF --indep-pairwise 50 10 $LD_parameter --make-bed --out $file\_LD
```

**Prune samples based on previous step**
```
plink --bfile $file\_HWE_MAF --extract $file\_LD.prune.in --make-bed --out $file\_HWE_MAF_LD
echo 'Pruned samples based on LD'
```

&nbsp;

## Generating PCA
**Make PCA and report number of SNPs used**
```
plink --bfile $file\_HWE_MAF_LD --pca 4 --out $file\_PCA
echo 'Generated eigenvalues for PCA'
SNP_count=$(cat BeanGoose_PCA.log | grep "pass" | awk '{print $1;}')
echo "PCA based on $SNP_count SNPs"
```

&nbsp;

## Plotting in R
Based on the R-script Plot_PCA.R
```R
#!/usr/bin/env

# Libraries
library(ggplot2)

####################
### Loading data ###
####################

# Load and process data
pca <- read.table("BeanGoose_PCA.eigenvec") # Load data
colnames(pca) <- c("Individual","X","PC1", "PC2", "PC3", "PC4") # Change column names

# Add column with species
pca$Species <- c(rep("A_Fabalis", 9),
                 rep("A_Serrirostris",9))
pca$Species <- factor(pca$Species)

################
### Plotting ###
################

# Choose colors
my_cols <- c("red", "blue")

# Set species names
species_names <- c("Taiga Bean Goose (A. fabalis)",
                     "Tundra Bean Goose (A. serrirostris)")

# Open PDF to save plot
pdf("PCA.pdf")

# Plot PCA
ggplot(data = pca, aes(x=PC1, y=PC2, col=Species)) +
  geom_point(size = 3) +
  theme_bw() +
  scale_color_manual(
    values = my_cols,
    labels = species_names,
    name = "") +
  theme(legend.position = "bottom")

# Close the pdf file
dev.off()
```

# Window Analyses
I will calculate different summary statistics (Fst, pi and Dxy) for a variety of window sizes (10kb, 20kb, 50kb, 100kb and 200kb). 
The distribution of these statistics will be plotted with R and I will check for correlations between them.

## Calculating Fst
To estimate Fst, I use VCFtools version 0.1.15.
The line of code below only takes bi-allelic SNPs and removes any variants with missing data.
```
vcftools --vcf BeanGoose.SNPs.HardFilt.vcf.gz --weir-fst-pop Fabalis.txt --weir-fst-pop Rossicus.txt --fst-window-size [X] --fst-window-step [X] --max-missing-count 0 --max-alleles 2 --min-alleles 2 —remove-filtered-all
```

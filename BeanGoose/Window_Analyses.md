# Window Analyses
I will calculate different summary statistics (Fst, pi and Dxy) for a variety of window sizes (10kb, 20kb, 50kb, 100kb and 200kb). 
The distribution of these statistics will be plotted with R and I will check for correlations between them.

&nbsp;

## Calculating Fst
To estimate Fst, I use VCFtools version 0.1.15.
The line of code below only takes bi-allelic SNPs and removes any variants with missing data.
```
vcftools --vcf BeanGoose.SNPs.HardFilt.vcf.gz --weir-fst-pop Fabalis.txt --weir-fst-pop Rossicus.txt --fst-window-size [X] --fst-window-step [X] --max-missing-count 0 --max-alleles 2 --min-alleles 2 —remove-filtered-all
```

&nbsp;

## Calculating pi and Dxy
To estimate pi and Dxy, I use the script by Simon Martin which can be found [here](https://github.com/simonhmartin/genomics_general).
First, I convert the VCF-file into the proper format.
```
python parseVCF.py -i BeanGoose.SNPs.HardFilt.vcf --skipIndels --minQual 30 --gtf flag=DP min=10 | gzip > BeanGoose.geno.gz
```
Next, I can calculate the summary statistics.
```
python popgenWindows.py --windType coordinate -w 10000 -m 10 -g BeanGoose.geno.gz -o BeanGoose_10kb.csv -T 8 -f phased -p Fabalis -p Rossicus --popsFile Populations.txt
```

&nbsp;

## Plotting in R
Here is the R-code to make violin-plots of the statistics.

And the R-code to plot and check for correlations.

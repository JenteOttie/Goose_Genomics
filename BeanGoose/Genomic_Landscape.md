# Genomic Landscape

I visualized the distribution of three summary statistics (Fst, Dxy and Pi) across the genome. 
I opted for 200kb windows, because smaller window sizes make the resulting picture to chaotic.
Because there is no chromosome-level assembly available for geese, I mapped all scaffolds to the Chicken (*Gallus gallus*) genome.
All chromosomes are thus pseudo-chromosomes.

&nbsp;

## Calculating summary statistics
This was done in the previous step of this project, [Window Analyses](https://github.com/JenteOttie/Goose_Genomics/edit/master/BeanGoose/Window_Analyses.md)
But for clarity, I will repeat the procedure here.

To estimate Fst, I use VCFtools version 0.1.15. As mentioned above, window size is 200kb
The line of code below only takes bi-allelic SNPs and removes any variants with missing data. 
```
vcftools --vcf BeanGoose.SNPs.HardFilt.vcf --weir-fst-pop Fabalis.txt --weir-fst-pop Rossicus.txt --fst-window-size 200000 --fst-window-step 200000 --max-missing-count 0 --max-alleles 2 --min-alleles 2 —remove-filtered-all
```

To estimate pi and Dxy, I use python-scripts written by Simon Martin, which can be found [here](https://github.com/simonhmartin/genomics_general). Important to keep in mind: these scripts need Python2.7 (running it with newer versions does not work)

First, I convert the VCF-file into the proper format, only using SNPs with a minumum quality of 30 and a read depth of 10x.
```
python parseVCF.py -i BeanGoose.SNPs.HardFilt.vcf --skipIndels --minQual 30 --gtf flag=DP min=10 | gzip > BeanGoose.geno.gz
```
Next, I can calculate the summary statistic for 200kb windows with a minimum of 100 SNPs per window.
```
python popgenWindows.py --windType coordinate -w 200000 -m 100 -g BeanGoose.geno.gz -o BeanGoose_200kb.csv -T 8 -f phased -p Fabalis -p Rossicus --popsFile Populations.txt
```

&nbsp;

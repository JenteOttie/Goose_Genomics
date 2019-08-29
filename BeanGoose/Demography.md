# Demographic Analyses
To infer the evolutionary history of Taiga and Tundra Bean Goose, I will do demographic modelling using [dadi](https://bitbucket.org/gutenkunstlab/dadi/src/master/).
This software models allele frequency changes using diffusion approximation. For more information about this approach, you can read this [paper](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000695).

![Dadi](https://journals.plos.org/plosgenetics/article/figure/image?size=large&id=10.1371/journal.pgen.1000695.g001)

## Preparing data
Because demographic modelling assumes neutral evolution, it is advised to use non-coding SNPs. I will select these SNPs using the software [snpEff](http://snpeff.sourceforge.net/).

**Installing snpEff**
```
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
```

**Create database for Swan Goose genome**
The software snpEff has several databases installed. You can check them with the command.
```
$ java -jar snpEff.jar databases
```
Because the Swan Goose (Anser cygnoides) genome is not installed, I need to create it myself. First, I adjust the file snpEff.config. I add the following lines in the database section:
```
# Swan Goose genome, version 1
ansCyg1.genome : Swan Goose
```
Then I make folders data/ and genome/, which will hold the data from the Swan Goose genome. In the data-folder I save the gff-file and in the genome-folder the fasta-file. I change the names of these files to genes.gff.gz and ansCyg1.fa.gz.
Download data from https://www.ncbi.nlm.nih.gov/genome/?term=txid8845[orgn]
```
$ mv GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.gff.gz genes.gff.gz
$ mv GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.fna.gz ansCyg1.fa.gz
```
Finally, I can build the database for this genome
```
$ java -jar snpEff.jar build -gff3 -v ansCyg1
```

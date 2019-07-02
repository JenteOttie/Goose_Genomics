# Pipeline to Reconstruct the Goose Phylogeny

## 1. Make fasta-file for each sample

This step is done with the bash-script _Generate_Fasta.sh_ which takes a BAM-file as input. 
The BAM-file is recalibrated (BQSR), realigned (GATK) and duplicated-marked (Picard).

First a fastq-file is generated with samtools mpileup. The reading depth is set between 5X and 100X.

$ samtools mpileup -uf $fasta $1 | bcftools call -c | vcfutils.pl vcf2fq -d 5 -D 100 >$ind.fq

Then the fastq-file is converted to fasta (one-liner based on https://www.biostars.org/p/85929/)

$ cat $ind.fq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' >$ind.fa

&nbsp;

## 2. Select specific window for analysis

From the generated fasta-files, I select particular windows for the phylogenetic analyses. 
These windows are saved to the file "Raw_Sequences.fa" using the python-script _Make_Fasta_for_Phylogeny.py_.
This script randomly picks a window of 10kb and checks the quality by counting the number of Ns in this window. 
If the number of Ns is above a certain threshold (e.g., 25%), the window is discarded.

&nbsp;

## 3. Align sequences

The sequences are aligned using Muscle version 3.8.31 (part of bash-script _Phylogeny_Maker.sh_)

$ muscle -in Raw_Sequences.fa -out Alignment.fa

&nbsp;

## 4. Model selection and phylogenetic analysis

The alignment "Alignment.fa" is used for the phylogenetic analysis with IQtree version 1.6.10 (part of bash-script _Phylogeny_Maker.sh_)
First, a model selection is run (option -m), then a phylogeny is estimated with 1000 ultrafast bootstraps (option -bb).
The -nt option chooses the best number of cores to run in parallel.

$ iqtree -s Alignment.fa -m MFP -bb 1000 -nt AUTO

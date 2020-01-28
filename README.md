# GooseGenomics
This repository contains the workflow and scripts that I used throughout my goose work at Uppsala University.

## Sample Preparation
I collected blood and tissue samples from several goose species. Genomic DNA was isolated using the Qiagen Gentra kit (Qiagen Inc.). Quality and quantity of the DNA was measured using the Qubit (Invitrogen, Life Technologies). Sequencing libraries were prepared from 100 ng DNA using the TruSeq Nano DNA sample preparation kit (cat# FC-121-4001/4002, Illumina Inc.), targeting an insert size of 350 bp. Paired-end sequencing (150 bp) was performed on an Illumina HiSeqX following standard procedures. 

Here is an overview of the samples
- Greater White-fronted Goose (*Anser albifrons*) - 10 samples
- Greylag Goose (*Anser anser*) - 13 samples
- Pink-footed Goose (*Anser brachyrhynchus*) - 15 samples
- Lesser White-fronted Goose (*Anser erythropus*) - 3 samples
- Taiga Bean Goose (*Anser fabalis*) - 9 samples
- Tundra Bean Goose (*Anser serrirostris*) - 9 samples
- Brent Goose (*Branta bernicla*) - 5 samples
- Canada Goose (*Branta canadensis*) - 2 samples
- Barnacle Goose (*Brana leucopsis*) - 5 samples

## Genomic Analyses
These samples were all mapped to the Swan Goose (*Anser cygnoides*) genome and SNPs were called using the GATK best practises. Resulting BAM-files and GVCF-files have been stored on the UPPMAX supercomputer. This (page)[https://github.com/JenteOttie/Goose_Genomics/blob/master/Data_Processing.md] walks you through the whole process

This project involved two main topics:
- Evolution of Taiga and Tundra Bean Goose (see here)[https://github.com/JenteOttie/Goose_Genomics/tree/master/BeanGoose]
- Phylogenomics of the genus *Anser* (see here)[https://github.com/JenteOttie/Goose_Genomics/tree/master/Phylogenomics[

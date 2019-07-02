#!/bin/bash -l
#SBATCH -A p2018002
#SBATCH -J fasta_generator
#SBATCH -p node
#SBATCH -n 8
#SBATCH -t 20:00:00

# Script to convert BAM-file into fasta-file
# Written by Jente Ottenburghs on 2019-01-16

# load modules
module load bioinfo-tools
module load samtools/1.8
module load bcftools/1.8

#Input files
filename=$1
ind=$(basename $filename .recal.realn.marked.bam)
fasta=/proj/sllstore2017033/nobackup/work/jente/Reference_Genome/ansCyg.fa

# Make fastq file
echo 'Making fastq-file'
date
samtools mpileup -uf $fasta $1 | bcftools call -c | vcfutils.pl vcf2fq -d 5 -D 100 >$ind.fq
echo 'Finished fastq-file'
date

# Convert fastq to fasta (one-liner based on https://www.biostars.org/p/85929/)
cat $ind.fq | paste - - - - | sed 's/^@/>/g'| cut -f1-2 | tr '\t' '\n' >$ind.fa

# Remove fastq-file
rm $ind.fq

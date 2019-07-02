#!/bin/bash -l
#SBATCH -A p2018002
#SBATCH -J Phylogeny
#SBATCH -p node
#SBATCH -n 8
#SBATCH -t 2:00:00

# Script to make phylogenetic trees from genome-files
# Written by Jente Ottenburghs on 2019-07-02

# Workflow
# 1. Extract windows from genomes (see Make_Fasta_for_Phylogeny.py)
# 2. Align windows with Muscle
# 3. Determine evolutionary model with IQ-tree
# 4. Make phylogenetic tree

# Make input-file (fasta)
python Make_Fasta_for_Phylogeny.py
echo 'Made input-file (fasta) for phylogenetic analyses)'
date

# Align sequences
ml bioinfo-tools
ml muscle/3.8.31 # Load module
muscle -in Raw_Sequences.fa -out Alignment.fa
echo 'Finished alignment with Muscle'
date

# Determine evolutionary model with IQ-tree
ml iqtree/1.6.10-omp-mpi # Load module
iqtree -s Alignment.fa -m MFP -bb 1000 -nt AUTO # Run IQ-tree with model selection (-m) and 1000 ultrafast bootstraps (-bb). The -nt option chooses the best number of cores to run in parallel.

# Delete all files, expect output-tree
cat Alignment.fa.treefile >>GeneTrees.txt
rm Alignment.f*

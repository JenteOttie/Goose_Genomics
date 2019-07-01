# From raw Fastq-files to variants in VCF-files
## Mapping reads to reference genome
All reads will be mapped to the Swan Goose (*Anser cygnoides*) genome.

Location on Uppmax: /proj/sllstore2017033/nobackup/work/jente/Reference_Genome/ansCyg.fa

More information about this genome at NCBI: https://www.ncbi.nlm.nih.gov/assembly/GCF_000971095.1

### Workflow
Make index for reference genome

$ samtools faidx ansCyg.fa

Make temporary directory to store files

$ mkdir

Go to sample directory (= input when running bash-script)

$ cd $1


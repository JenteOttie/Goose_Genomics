# From raw Fastq-files to variants in VCF-files
## 1. Mapping reads to reference genome
All reads will be mapped to the Swan Goose (*Anser cygnoides*) genome.

Location on Uppmax: /proj/sllstore2017033/nobackup/work/jente/Reference_Genome/ansCyg.fa

More information about this genome at NCBI: https://www.ncbi.nlm.nih.gov/assembly/GCF_000971095.1

&nbsp;

**Make index for reference genome**

$ samtools faidx ansCyg.fa

**Go to sample directory (= input when running bash-script)**

$ cd $1

**Make temporary directory and copy fastq-files to workplace for analyses**

$ mkdir /proj/sllstore2017033/nobackup/work/jente/temp_$ind

$ cp \*.gz /proj/sllstore2017033/nobackup/work/jente/temp_$ind

**Go to work-directory and unzip files**

$cd /proj/sllstore2017033/nobackup/work/jente/temp_$ind

$ \gunzip *.gz

**Set file locations**

$ ref=/proj/sllstore2017033/nobackup/work/jente/Reference_Genome/ansCyg

$ fq1=\*R1_001.fastq

$ fq2=\*R2_001.fastq

**Run BWA (version 0.7.17)**

$ bwa mem -R "@RG\tID:$ind\tSM:$ind\tLB:$ind\tPI:350\tPL:Illumina" -t 20 -M $ref $fq1 $fq2 | samtools sort -T $SNIC_TMP/$ind - >$ind.bam.tmp && mv $ind.bam.tmp $ind.bam

**Calculate mapping statistics with samtools flagstat**

$ samtools flagstat $ind.bam >$ind.stats.txt

**Move files to BAM-file directory**

$ mv \*.bam /proj/sllstore2017033/nobackup/work/jente/BAM_Files

## 2. Marking duplicates with PICARD (version 2.10.3)

$ filename=$1

$ ind=$(basename $filename .bam)

$ java -jar $PICARD_HOME/picard.jar MarkDuplicates I=$1 O=$ind.marked.bam METRICS_FILE=$ind.marked.metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true

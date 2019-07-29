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

$ cd /proj/sllstore2017033/nobackup/work/jente/temp_$ind

$ \gunzip *.gz

**Set file locations**

$ ref=/proj/sllstore2017033/nobackup/work/jente/Reference_Genome/ansCyg

$ fq1=\*R1_001.fastq

$ fq2=\*R2_001.fastq

**Run BWA (version 0.7.17) and output BAM-file**

$ bwa mem -R "@RG\tID:$ind\tSM:$ind\tLB:$ind\tPI:350\tPL:Illumina" -t 20 -M $ref $fq1 $fq2 | samtools sort -T $SNIC_TMP/$ind - >$ind.bam.tmp && mv $ind.bam.tmp $ind.bam

**Calculate mapping statistics with samtools flagstat**

$ samtools flagstat $ind.bam >$ind.stats.txt

**Move files to BAM-file directory**

$ mv \*.bam /proj/sllstore2017033/nobackup/work/jente/BAM_Files

&nbsp;

## 2. Marking duplicates with PICARD (version 2.10.3)

$ filename=$1

$ ind=$(basename $filename .bam)

$ java -jar $PICARD_HOME/picard.jar MarkDuplicates I=$1 O=$ind.marked.bam METRICS_FILE=$ind.marked.metrics REMOVE_DUPLICATES=false ASSUME_SORTED=true

&nbsp;

## 3. Local Realignment with GATK (version 3.7)

**Set file locations and names**

$ filename=$1

$ ind=$(basename $filename .marked.bam)
$ fasta=/proj/sllstore2017033/nobackup/work/jente/Reference_Genome/ansCyg.fa

**Create index of Bam-file**

$ samtools index $1

**Make intervals**

$ java -Xmx60g -jar $GATK_HOME/GenomeAnalysisTK.jar -T RealignerTargetCreator -I $1 -R $fasta -o $ind.intervals


**Realignment (this step cannot run in parallel)**

$ java -Xmx60g -Djava.io.tmpdir=$SNIC_TMP -jar $GATK_HOME/GenomeAnalysisTK.jar -T IndelRealigner -I $1 -R $fasta -targetIntervals $ind.intervals -o $ind.realn.marked.bam.tmp && mv $ind.realn.marked.bam.tmp $ind.realn.marked.bam && mv $ind.realn.marked.bam.tmp.bai $ind.realn.marked.bam.bai

&nbsp;

## 4. Base Quality Score Recalibration (BQSR) with GATK (version 3.7)

Normally, you would use a set of high quality SNPs to do the BQSR. Because this is not available for geese, I will perform a bootstrapping approach. This strategy is explained in more detail at the [GATK website](https://software.broadinstitute.org/gatk/documentation/article?id=11081). This is done on an indidivual level.

**Set file locations and names**

$ BAM=$1

$ ind=$(basename $1 .realn.marked.bam)

$ fasta=/proj/sllstore2017033/nobackup/work/jente/Reference_Genome/ansCyg.fa

**Run GATK HaplotypeCaller**

$ java -Xmx120g -Djava.io.tmpdir=$SNIC_TMP -jar $GATK -T HaplotypeCaller -R $fasta -I $BAM -nct 20 -o $ind.raw.snps.indels.vcf

$ VCF=$ind.raw.snps.indels.vcf

**Extract SNPs with minimum depth of 10 (DP>10) and MQRankSum<-0.2 (based on script of Pall)**

$ java -Xmx120g -Djava.io.tmpdir=$SNIC_TMP -jar $GATK -T SelectVariants -R $fasta -V $VCF -selectType SNP -select "MQRankSum<-0.2 && DP>10" -o $ind.recal.vcf

**Make recalibration tables**

$ java -Xmx120g -Djava.io.tmpdir=$SNIC_TMP -jar $GATK -T BaseRecalibrator -R $fasta -I $BAM -knownSites $ind.recal.vcf -nct 20 -o $ind.recal1.table

$ java -Xmx120g -Djava.io.tmpdir=$SNIC_TMP -jar $GATK -T BaseRecalibrator -R $fasta -I $BAM -knownSites $ind.recal.vcf -BQSR $ind.recal1.table -nct 20 -o $ind.recal2.table

**Compare recalibration tables (output to PDF-file)**

$ java -jar $GATK -T AnalyzeCovariates -R $fasta -before $ind.recal1.table -after $ind.recal2.table -plots $ind.BQSR.1vs2.pdf

**Recalibrate BAM-file and move to BAM-folder**

$ java -Xmx120g -Djava.io.tmpdir=$SNIC_TMP -jar $GATK -T PrintReads -R $fasta -I $BAM -BQSR $ind.recal2.table -nct 20 -o BQRS_PDF/$ind.recal.realn.marked.bam

$ mv $ind.recal.realn.marked.b* BAM_Files/

&nbsp;

## 5. Joint genotyping with GATK (version 3.7)

All individual VCF-file are combined in one file using the function GenotypeGVCFs

$ java -Xmx7g -jar $GATK -T GenotypeGVCFs -R $fasta \
--variant [ind1] \
--variant [ind2] \
... \
-o [output.vcf]

&nbsp;

## 6. Apply hard-filter to VCF-file with GATK (version 3.7)

**Set file locations and names**

$ filename=$1

$ ind=$(basename $filename .vcf)

$ fasta=/proj/sllstore2017033/nobackup/work/jente/Reference_Genome/ansCyg.fa

**Extract SNPs from file**
java -Xmx6g -jar $GATK -T SelectVariants -R $fasta -V $1 -selectType SNP -o $ind.SNPs.vcf

**Run hard SNP filter on VCF-file - Setting based on GATK best practises**

$ java -Xmx6g -jar $GATK -T VariantFiltration -R $fasta -V $ind.SNPs.vcf --filterExpression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" --filterName "hard_filt" -o $ind.SNPs.HardFilt.vcf

&nbsp;

## 7. Downstream analyses

The final VCF-file can be used for further analyses. Different filtering steps will be applied depending on the analysis.

# Genomic Landscape

I visualized the distribution of three summary statistics (Fst, Dxy and Pi) across the genome. 
I opted for 200kb windows, because smaller window sizes make the resulting picture to chaotic.
Because there is no chromosome-level assembly available for geese, I mapped all scaffolds to the Chicken (*Gallus gallus*) genome.
All chromosomes are thus pseudo-chromosomes.

&nbsp;

## Calculating summary statistics
This was done in the previous step of this project, [Window Analyses](https://github.com/JenteOttie/Goose_Genomics/edit/master/BeanGoose/Window_Analyses.md).
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

## Mapping Scaffolds to Chicken Genome
I used LASTZ version 1.04.00 to map all goose scaffolds to the Chicken genome. This software takes two fasta-files as input and outputs the location of all scaffolds on the genome. the parameters were tuned by Linnea Smeds and Alexander Suh for aligning bird species.
```
GOOSE=ansCyg.fa
REF=galGal6.fa

lastz_32 $GOOSE[multiple] $REF M=254 K=4500 L=3000 Y=15000 C=2 T=2 --matchcount=10000 --ambiguous=iupac --format=general:name1,start1,end1,length1,strand1,name2,start2,end2,length2,strand2
```

&nbsp;

## Creating the Genomic Landscape
The files generated in the previous steps (summary statistics in 200kb windows and LASTZ-alignment) are combined in an R-script to create a genomic landscape.
```R
library(ggplot2) # plotting
library(Rmisc) # multiplot
library(cowplot) # For function plotgrid

#################
### Load data ###
#################

# output from LASTZ
mapping <- read.table("LASTZ_Chicken.output")
colnames(mapping) <- c("scaffold", "start1", "end1", "length1", "strand1", "chromosome", "start2", "end2", "length2", "strand2")

# Fst data
Fst_Fa_Ro <- read.table("BeanGoose_200kb.windowed.weir.fst", header = T )
colnames(Fst_Fa_Ro) <- c("scaffold", "start", "end", "N_variants", "Weighted_Fst", "Mean_Fst")
Fst_Fa_Ro$Mean_Fst[Fst_Fa_Ro$Mean_Fst < 0] <- 0

# Dxy and Pi data
Diversity <- read.csv("BeanGoose_200kb.csv", header = TRUE)
Dxy_Fa_Ro <- Diversity[,c(1,2,3,10)]
colnames(Dxy_Fa_Ro) <- c("scaffold", "start", "end", "Dxy")

Pi_Fa <- Diversity[,c(1,2,3,6)]
colnames(Pi_Fa) <- c("scaffold", "start", "end", "Pi")

Pi_Ro <- Diversity[,c(1,2,3,7)]
colnames(Pi_Ro) <- c("scaffold", "start", "end", "Pi")

# Chicken chromosomes (contains only chromosomes with sufficient scaffolds mapped to it)
chicken_chr_subset <- read.table("Chicken_chr.txt", header = T )

chromosome_numbers <- chicken_chr$chr

for (c in 1:length(chicken_chr$chr)){

  #######################
  ### Order scaffolds ###
  #######################

  # Select chromosome
  chr <- chicken_chr$chr[c]

  # Extract data for chromosome
  chromosome = as.character(chicken_chr[which(chicken_chr$chr == chr),][,2]) # Get RefSeq-name for chromosome
  chr <- mapping[which(mapping$chromosome == chromosome),]

  # Order starting positions on chromosome (start2)
  chr_ordered <- chr[order(chr$start2),]

  # Extract scaffolds
  # The scaffolds are now ordered according to their starting positions
  scaffolds <- unique(chr_ordered$scaffold)

  ################
  ### Plot Fst ###
  ################

  # This loop takes the Fst-data per scaffold and orders the scaffolds according to the data above.
  Fst_Fa_Ro_chr <- {}
  Extra_bp = 0
  for (i in 1:length(scaffolds)){

    # Get scaffold name
    scaffold_name = as.character(scaffolds[i])
    scaffold_Fst_Fa_Ro <- Fst_Fa_Ro[which(Fst_Fa_Ro$scaffold == scaffold_name),]
    
    # Get value of final window to update the bp later
    # I need to get this value before adding the extra bp. Otherwise it will increase exponentially
    final_end <- tail(scaffold_Fst_Fa_Ro$end, 1)
  
    # Add extra bp to start and end
    scaffold_Fst_Fa_Ro$start <- scaffold_Fst_Fa_Ro$start + Extra_bp
    scaffold_Fst_Fa_Ro$end <- scaffold_Fst_Fa_Ro$end + Extra_bp

    # Add data to dataframe
    Fst_Fa_Ro_chr <-rbind(Fst_Fa_Ro_chr, scaffold_Fst_Fa_Ro)
  
    # Update extra bp
    Extra_bp = Extra_bp + final_end
  }

  # Plot genomic landscape
  p_Fst <- ggplot(Fst_Fa_Ro_chr, aes(x = start/1000, y = Mean_Fst)) +
    geom_area(fill = 'lightblue', color = 'black') +
    theme_bw() +
    ylab("Fst") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_text(size=7),
          axis.text.y=element_text(size=5)) +
    ggtitle(chromosome_numbers[c]) +
    scale_y_continuous(limits = c(0, 0.5),
                       breaks = c(0, 0.25, 0.5))

    
  ################
  ### Plot Dxy ###
  ################
  
  # This loop takes the Dxy-data per scaffold and orders the scaffolds according to the data above.
  Dxy_Fa_Ro_chr <- {}
  Extra_bp = 0
  for (i in 1:length(scaffolds)){
    
    # Get scaffold name
    scaffold_name = as.character(scaffolds[i])
    scaffold_Dxy_Fa_Ro <- Dxy_Fa_Ro[which(Dxy_Fa_Ro$scaffold == scaffold_name),]
    
    # Get value of final window to update the bp later
    # I need to get this value before adding the extra bp. Otherwise it will increase exponentially
    final_end <- tail(scaffold_Dxy_Fa_Ro$end, 1)
    
    # Add extra bp to start and end
    scaffold_Dxy_Fa_Ro$start <- scaffold_Dxy_Fa_Ro$start + Extra_bp
    scaffold_Dxy_Fa_Ro$end <- scaffold_Dxy_Fa_Ro$end + Extra_bp
    
    # Add data to dataframe
    Dxy_Fa_Ro_chr <- rbind(Dxy_Fa_Ro_chr, scaffold_Dxy_Fa_Ro)
    
    # Update extra bp
    Extra_bp = Extra_bp + final_end
  }
  
  # Plot genomic landscape
  p_Dxy <- ggplot(Dxy_Fa_Ro_chr, aes(x = start/1000, y = Dxy)) +
    geom_area(fill = 'lightblue', color = 'black') +
    theme_bw() +
    ylab("Dxy") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_text(size=7),
          axis.text.y=element_text(size=5)) +
    scale_y_continuous(limits = c(0, 0.5),
                       breaks = c(0, 0.25, 0.5))

  ###############
  ### Plot Pi ###
  ###############
  
  # This loop takes the Pi-data per scaffold and orders the scaffolds according to the data above.
  # This analysis has to be done separately for A. fabalis and A. serrirostris
  
  # Fabalis
  Pi_Fa_chr <- {}
  Extra_bp = 0
  for (i in 1:length(scaffolds)){
    
    # Get scaffold name
    scaffold_name = as.character(scaffolds[i])
    scaffold_Pi_Fa <- Pi_Fa[which(Pi_Fa$scaffold == scaffold_name),]
    
    # Get value of final window to update the bp later
    # I need to get this value before adding the extra bp. Otherwise it will increase exponentially
    final_end <- tail(scaffold_Pi_Fa$end, 1)
    
    # Add extra bp to start and end
    scaffold_Pi_Fa$start <- scaffold_Pi_Fa$start + Extra_bp
    scaffold_Pi_Fa$end <- scaffold_Pi_Fa$end + Extra_bp
    
    # Add data to dataframe
    Pi_Fa_chr <- rbind(Pi_Fa_chr, scaffold_Pi_Fa)
    
    # Update extra bp
    Extra_bp = Extra_bp + final_end
  }
  
  # Rossicus
  Pi_Ro_chr <- {}
  Extra_bp = 0
  for (i in 1:length(scaffolds)){
    
    # Get scaffold name
    scaffold_name = as.character(scaffolds[i])
    scaffold_Pi_Ro <- Pi_Ro[which(Pi_Ro$scaffold == scaffold_name),]
    
    # Get value of final window to update the bp later
    # I need to get this value before adding the extra bp. Otherwise it will increase exponentially
    final_end <- tail(scaffold_Pi_Ro$end, 1)
    
    # Add extra bp to start and end
    scaffold_Pi_Ro$start <- scaffold_Pi_Ro$start + Extra_bp
    scaffold_Pi_Ro$end <- scaffold_Pi_Ro$end + Extra_bp
    
    # Add data to dataframe
    Pi_Ro_chr <- rbind(Pi_Ro_chr, scaffold_Pi_Ro)
    
    # Update extra bp
    Extra_bp = Extra_bp + final_end
  }
  
  # Put the dataframes of Pi_Fa and Pi_Ro together and add column with species names
  Pi_Fa_Ro_chr <- rbind(Pi_Fa_chr, Pi_Ro_chr)
  Pi_Fa_Ro_chr$Species <- c(rep("Fabalis", length(Pi_Ro_chr$scaffold)), rep("Rossicus", length(Pi_Fa_chr$scaffold)))
  
  # Plot genomic landscape
  p_Pi <- ggplot(Pi_Fa_Ro_chr, aes(x = start/1000, y = Pi, color = Species)) +
    geom_line() +
    theme_bw() +
    ylab(expression(pi)) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_text(size=7),
          axis.text.y=element_text(size=5)) +
    theme(legend.position="none") +
    scale_y_continuous(limits = c(0, 0.5),
                       breaks = c(0, 0.25, 0.5))
  
  ###########################
  ### Put it all together ###
  ###########################
  
  p <- plot_grid(p_Fst, p_Dxy, p_Pi, nrow=3, rel_heights = c(1.5,1,1))
  
  # Give plot a name that corresponds to chromosome number
  plot_name <- paste("p", chromosome_numbers[c], sep = "")
  assign(plot_name,p)
}

layout(matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1,
                2, 2, 2, 2, 2, 2, 2, 2, 2,
                3, 3, 3, 3, 3, 4, 4, 4, 4,
                5, 5, 5, 6, 6, 7, 7, 8, 8,
                9, 9, 10, 10, 11, 12, 13, 14, 15,
                16, 17, 18, 19, 20, 21, 22, 23, 24,
                25, 26, 27, 27, 27, 28, 28, 28, 28
                ), nrow=7, byrow=TRUE))
layout.show(n=28)


multiplot(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10,
          p11, p12, p13, p14, p15, p17, p18, p19, p20,
          p21, p22, p23, p24, p26, p27, p28, pW, pZ,
          layout=matrix(c(1, 1, 1, 1, 1, 1, 1, 1, 1,
                   2, 2, 2, 2, 2, 2, 2, 2, 2,
                   3, 3, 3, 3, 3, 4, 4, 4, 4,
                   5, 5, 5, 6, 6, 7, 7, 8, 8,
                   9, 9, 10, 10, 11, 12, 13, 14, 15,
                   16, 17, 18, 19, 20, 21, 22, 23, 24,
                   25, 26, 27, 27, 28, 28, 28, 28, 28
          ), nrow=7, byrow=TRUE))
```

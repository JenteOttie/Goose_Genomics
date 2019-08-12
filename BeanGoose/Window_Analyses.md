# Window Analyses
I will calculate different summary statistics (Fst, pi and Dxy) for a variety of window sizes (10kb, 20kb, 50kb, 100kb and 200kb). 
The distribution of these statistics will be plotted with R and I will check for correlations between them.

&nbsp;

## Calculating Fst
To estimate Fst, I use VCFtools version 0.1.15.
The line of code below only takes bi-allelic SNPs and removes any variants with missing data.
```
vcftools --vcf BeanGoose.SNPs.HardFilt.vcf --weir-fst-pop Fabalis.txt --weir-fst-pop Rossicus.txt --fst-window-size [X] --fst-window-step [X] --max-missing-count 0 --max-alleles 2 --min-alleles 2 —remove-filtered-all
```

&nbsp;

## Calculating pi and Dxy
To estimate pi and Dxy, I use python-scripts written by Simon Martin, which can be found [here](https://github.com/simonhmartin/genomics_general). Important to keep in mind: these scripts need Python2.7 (running it with newer versions does not work)

First, I convert the VCF-file into the proper format.
```
python parseVCF.py -i BeanGoose.SNPs.HardFilt.vcf --skipIndels --minQual 30 --gtf flag=DP min=10 | gzip > BeanGoose.geno.gz
```
Next, I can calculate the summary statistics. Example here is for 10kb windows with a minimum of 100 SNPs per window.
```
python popgenWindows.py --windType coordinate -w 10000 -m 100 -g BeanGoose.geno.gz -o BeanGoose_10kb.csv -T 8 -f phased -p Fabalis -p Rossicus --popsFile Populations.txt
```

&nbsp;

## Plotting in R
Here is the R-code to make violin-plots and test for correlations between the summary statistics (for 200 kb windows). 
```R
library(ggplot2) # plotting
library(Rmisc) # multiplot

#################
### LOAD DATA ###
#################

# Load Fst-results
Fst <- read.table("G:/Data/Postdoc Uppsala/Bean Goose Project/Window Analyses/BeanGoose_200kb.windowed.weir.fst", header = T )
Fst <- Fst[,c(1,2,3,6)]
colnames(Fst) <- c("scaffold", "start", "end", "Fst")
Fst$Fst[Fst$Fst < 0] <- 0 # Turn negative Fst-values into zeros

Pi_Dxy <- read.csv("G:/Data/Postdoc Uppsala/Bean Goose Project/Window Analyses/BeanGoose_200kb.csv", header = TRUE)
Pi_Dxy <- Pi_Dxy[,c(1,2,3,6,7,8)]
colnames(Pi_Dxy) <- c("scaffold", "start", "end", "Pi_Fa", "Pi_Ro", "Dxy")

All_Data <- merge(Fst, Pi_Dxy, by=c("scaffold", "start", "end"))
All_Data$Extra <- rep(1, length(All_Data$Fst)) # Add extra column to make violin-plots

####################
### VIOLIN PLOTS ###
####################

# Fst
Fst <- ggplot(All_Data, aes(x=Extra, y=Fst)) +
  geom_violin(fill="lightblue") +
  geom_boxplot(width = 0.05) +
  theme_bw() +
  ylab("Fst") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Dxy
Dxy <- ggplot(All_Data, aes(x=Extra, y=Dxy)) +
  geom_violin(fill="lightblue") +
  geom_boxplot(width = 0.05) +
  theme_bw() +
  ylab("Dxy") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Pi (Anser fabalis)
Pi_Fa <- ggplot(All_Data, aes(x=Extra, y=Pi_Fa)) +
  geom_violin(fill="lightblue") +
  geom_boxplot(width = 0.05) +
  theme_bw() +
  ylab("\u03C0 (Anser fabalis)") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

# Pi (Anser rossicus)
Pi_Ro <- ggplot(All_Data, aes(x=Extra, y=Pi_Ro)) +
  geom_violin(fill="lightblue") +
  geom_boxplot(width = 0.05) +
  theme_bw() +
  ylab("\u03C0 (Anser serrirostris)") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

multiplot(Fst, Dxy, Pi_Fa, Pi_Ro, cols = 4)

####################
### CORRELATIONS ###
####################

# Fst vs. Dxy
correlation<-cor.test(All_Data$Fst, All_Data$Dxy, method = "spearman", exact = FALSE)
Fst_vs_Dxy <- ggplot(All_Data, aes(x=Fst, y=Dxy)) +
  geom_point() +
  theme_bw() +
  annotate('text', x=0.35, y=0.05, cex=2, label=paste("p-value = ", round(correlation$p.value, 27))) +
  annotate('text', x=0.35, y=0.03, cex=2, label=paste("\u03C1 = ", round(correlation$estimate, 2)))

# Fst vs. Pi_Fa
correlation<-cor.test(All_Data$Fst, All_Data$Pi_Fa, method = "spearman", exact = FALSE)
Fst_vs_Pi_Fa <- ggplot(All_Data, aes(x=Fst, y=Pi_Fa)) +
  geom_point() +
  theme_bw() +
  ylab("\u03C0 (Anser fabalis)") +
  annotate('text', x=0.35, y=0.05, cex=2, label=paste("p-value = ", round(correlation$p.value, 51))) +
  annotate('text', x=0.35, y=0.03, cex=2, label=paste("\u03C1 = ", round(correlation$estimate, 2)))

# Fst vs. Pi_Ro
correlation<-cor.test(All_Data$Fst, All_Data$Pi_Ro, method = "spearman", exact = FALSE)
Fst_vs_Pi_Ro <- ggplot(All_Data, aes(x=Fst, y=Pi_Ro)) +
  geom_point() +
  theme_bw() +
  ylab("\u03C0 (Anser serrirostris)") +
  annotate('text', x=0.35, y=0.05, cex=2, label=paste("p-value = ", round(correlation$p.value, 2))) +
  annotate('text', x=0.35, y=0.03, cex=2, label=paste("\u03C1 = ", round(correlation$estimate, 2)))

# Dxy vs. Pi_Fa
correlation<-cor.test(All_Data$Dxy, All_Data$Pi_Fa, method = "spearman", exact = FALSE)
Dxy_vs_Pi_Fa <- ggplot(All_Data, aes(x=Dxy, y=Pi_Fa)) +
  geom_point() +
  theme_bw() +
  ylab("\u03C0 (Anser fabalis)") +
  annotate('text', x=0.35, y=0.05, cex=2, label=paste("p-value = ", round(correlation$p.value, 2))) +
  annotate('text', x=0.35, y=0.03, cex=2, label=paste("\u03C1 = ", round(correlation$estimate, 2)))

# Dxy vs. Pi_Ro
correlation<-cor.test(All_Data$Dxy, All_Data$Pi_Ro, method = "spearman", exact = FALSE)
Dxy_vs_Pi_Ro <- ggplot(All_Data, aes(x=Dxy, y=Pi_Ro)) +
  geom_point() +
  theme_bw() +
  ylab("\u03C0 (Anser serrirostris)") +
  annotate('text', x=0.35, y=0.05, cex=2, label=paste("p-value = ", round(correlation$p.value, 2))) +
  annotate('text', x=0.35, y=0.03, cex=2, label=paste("\u03C1 = ", round(correlation$estimate, 2)))

# Pi_Fa vs. Pi_Ro
correlation<-cor.test(All_Data$Pi_Fa, All_Data$Pi_Ro, method = "spearman", exact = FALSE)
Pi_Fa_vs_Pi_Ro <- ggplot(All_Data, aes(x=Pi_Fa, y=Pi_Ro)) +
  geom_point() +
  theme_bw() +
  xlab("\u03C0 (Anser fabalis)") +
  ylab("\u03C0 (Anser serrirostris)") +
  annotate('text', x=0.35, y=0.05, cex=2, label=paste("p-value = ", round(correlation$p.value, 2))) +
  annotate('text', x=0.35, y=0.03, cex=2, label=paste("\u03C1 = ", round(correlation$estimate, 2)))

# Put it all together
multiplot(Fst_vs_Dxy, Dxy_vs_Pi_Fa, Fst_vs_Pi_Fa, Dxy_vs_Pi_Ro, Fst_vs_Pi_Ro, Pi_Fa_vs_Pi_Ro,
          cols = 3)
```

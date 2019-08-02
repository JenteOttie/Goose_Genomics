# Finding introgressed regions using a Machine Learning approach
## Background
In the previous step - demographic analysis - I showed that Taiga and Tundra Bean Goose have exchanged DNA during their evolutionary history.
To explore which genomic regions have been exchanged, I will apply a machine learning approach: FILET (Finding Introgressed Loci via Extra-Trees).
This software combines information from a suite of population genetic summary statistics that capture the genetic variation between two species (Schrider et al., 2018).
The machine learning algorithm relies on an Extra-Trees classifier (Geurts, Ernst, & Wehenkel, 2006), an extension of the Random Forest approach (Breiman, 2001). 
Random Forest is a collection of learning techniques that creates a forest of randomly generated decision trees that attempt to classify data points into predetermined groups. 
In this case, these groups correspond to introgression from population 1 into population 2, introgression from population 2 into population 1, and no introgression.

Here is a link to the FILET-paper: https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007341

&nbsp;

## Simulating training datasets
The algorithm needs to be trained. I will use simulated data as a training set. The data is simulated using [msmove](https://github.com/geneva/msmove).
Three scenarios are simulated: 
1. Introgression from population 1 to 2 
2. Introgression from population 2 to 1 
3. No Introgression
```
msmove 18 10000 -t 50 -I 2 9 9 -n 1 1.40176528 -n 2 0.6912292 -ma x 5.80351132 26.15397712 x -ev 0.14821583 1 2 1 >mig12.txt
msmove 18 10000 -t 50 -I 2 9 9 -n 1 1.40176528 -n 2 0.6912292 -ma x 5.80351132 26.15397712 x -ev 0.14821583 2 1 1 >mig21.txt
msmove 18 10000 -t 0.1 -I 2 9 9 0.0001 -n 1 1.40176528 -n 2 0.6912292 -ev 0.14821583 1 2 0 >Nomig.txt
```
The parameters are listed below. All numbers are based on the demographic model in the previous step.
- -t (theta)
- -I (number of populations) (samples in pop1) (samples in pop2)
- -n 1 (size of pop1)
- -n 2 (size of pop2)
- -ma (matrix with migration rates)
- -ev (time of migration) (source pop) (sink pop) (likelihood of migration event)

&nbsp;

## Preparing data for analyses
The input data for FILET consists of fasta-files for all populations and a bed-file with the coordinates of the windows to be analyzed.
First, I create the bed-file. This is done using a csv-file that I generated in the window-analyses (using the scripts by Simon Martin).
```
cut -d "," -f 1,2,3 10kb.csv >Test.txt      # Extract first three lines
sed -i '1d' Test.txt 		                   # Remove first line 
sed -i -e 's/,/\t/g' Test.txt		           # Replace comma (,) with tab (\t)
```
Next, I make fasta- and bed-files for each scaffold using the python-script [Prepare_Fasta_for_ML.py](https://github.com/JenteOttie/Goose_Genomics/blob/master/BeanGoose/Prepare_Fasta_for_ML.py)

&nbsp;

## Training the machine learning algoritm
Before I can analyse the data, I will need to train the algorithm with the simulated training sets. This is done using the python-scripts available in FILET. Because my server (UPPMAX) does not support some of the python packages, I first start a virtual environment where I can load and update my own python-packages.
```
source /home/jente/dadi/bin/activate
```
**Calculate summary statistics for simulated data (10kb windows)**
```
for inFile in `ls trainingSims/ | grep .msOut` ; do cat trainingSims/$inFile | ../FILET-master/twoPopnStats_forML 9 9 | python ../FILET-master/normalizeTwoPopnStats.py None 10000 >trainingSimsStats/$inFile; done
```
**Train classifier algorithm**
```
python buildThreeClassTrainingSet.py trainingSimsStats/ trainingSets/threeClass.fvec
python trainFiletClassifier.py trainingSets/threeClass.fvec classifier/threeClass.p pi1 hetVar1 ss1 private1 tajd1 ZnS1 pi2 hetVar2 ss2 private2 tajd2 ZnS2 Fst dxy_mean dxy_min gmin zx dd1 dd2 ddRank1 ddRank2
```

&nbsp;

## Classifying the data
Now it is time to classify the data. I do this on a scaffold-level. Each scaffold has its own fasta- and bed-files.
The results for each scaffold - summary stats and classification of windows - are saved in their own directory.
```
scaffold=$1                                   # Scaffold to classify
mkdir featureVectorsToClassify/${scaffold}/   # Make directory to save summary stats
mkdir results/${scaffold}/                    # Make directory to save results
```
**Calculate summary statistics on data**
```
/pgStatsBedSubpop_forML dataToClassify/${scaffold}_fabalis.fa dataToClassify/${scaffold}_rossicus.fa dataToClassify/${scaffold}_anc.fa dataToClassify/${scaffold}.bed 0.5 | python ../FILET-master/normalizePgStats.py 10000 >featureVectorsToClassify/${scaffold}/${scaffold}.ss
```
**Classify data**
```
python classifyChromosome.py classifier/threeClass.p featureVectorsToClassify/${scaffold}/ 0.05 results/${scaffold}/
```
All these commands can be found in the script Run_ML_pipeline.sh. To run it on all scaffolds in parallel, I use this line of code.
```
cat ../Scaffolds.txt | while read line; do sbatch Run_ML_pipeline.sh $line; done
```

&nbsp;

## Checking the results
To see which summary statistics were most important in classifying the windows, I run the python-script [getFeatureRankings.py](https://github.com/kern-lab/FILET/blob/master/getFeatureRankings.py) (written by Dan Schrider) on the file threeClass.p. This outputs a table that ranks the summary stats and shows their contribution to the classification. 

**Collecting the results**

The results for the different scaffolds have been saved in separate folders. To collect all the results in one file, just loop through the folders and extract the results.
```
for f in results/NW_01318*/* ; do cat $f >>Results.txt ; done
```
The results-file contains the following information (in separate columns)
1. Chromosome
2. Starting coordinate of the classified window
3. Ending coordinate
4. The number of sites in the window that are not 'N' across the entire sample (treating this as our number of informative sites)
5. The class label of the chosen class
6. The posterior probability of class zero
7. The posterior probability of class one
8. The posterior probability of class two

&nbsp;

**Analyzing the results: Gene content**

The number of windows per class can be calculated using this short Python-script (Sort_Introgressed_Windows.py)
```Python
import sys

# Open input and output-files
Windows_file = open(sys.argv[1], "r")
output_1 = open("Introgression_Class1.txt", "a")
output_2 = open("Introgression_Class2.txt", "a")

# Put introgressed windows in separate files
for line in Windows_file:
	# extract information about window and introgression class
	scaffold = line.split()[0]
	start = line.split()[1]
	end = line.split()[2]
	introgression_class = line.split()[4]
	
	window_coordinates = scaffold + "\t" + start + "\t" + end + "\n"
	
	if introgression_class == '1':
		output_1.write(window_coordinates)
		
	elif introgression_class == '2':
		output_2.write(window_coordinates)

# Close all files
Windows_file.close()
output_1.close()
output_2.close()
```
Next, I can check which genes are in the introgressed windows, using the script [Genes_in_Windows.py](https://github.com/JenteOttie/Goose_Genomics/blob/master/Genes_in_Windows.py). This script takes two inputs: the file with the windows and a GFF-file with the gene coordinates (which can be found [here](https://www.ncbi.nlm.nih.gov/genome/?term=txid8845[orgn]) for the goose genome). The resulting lists can be run through a [GO-term-analysis](http://geneontology.org/).

&nbsp;

**Analyzing the results: Summary Statistics**

To get an idea of the summary statistics for the different windows, we can plot them. First, we collect the summary stats data for all windows. The stats for all scaffolds can be found in the map featureVectorsToClassify. I put all the data in one file using this command.
```
for f in featureVectorsToClassify/NW_01318*/* ; do cat $f >>Window_Stats.txt ; done
```
Next, I use the Python-script Update_Introgressed_Windows.py to assign each window to the right introgression-class. This script takes two inputs: the file with all windows (Window_Stats.txt) and the file with the classified windows (Results.txt).
```Python
# Script to label introgressed regions based on Machine Learning analysis
# Written by Jente Ottenburghs on 2019-06-04

import sys

# Input files
# 1 = file with all windows
# 2 = file with introgressed windows

# Open input and output-file
Windows_file = open(sys.argv[1], "r")
output = open("Window_Stats_Classification.txt", "a")

# Put introgressed windows into a list. This will later be used to compare with the gene file.
Introgression_file = open(sys.argv[2], "r")
Introgression1 = []
Introgression2 = []
for line in Introgression_file:
	data = line.split()
	window = data[0] + ":" + data[1] + "-" + data[2]
	if data[4] == '1':
		Introgression1.append(window)
	elif data[4] == '2':
		Introgression2.append(window)
Introgression_file.close()

# Add header with all summary stats to file
output.write("chrom\tchromStart\tchromEnd\tnumSites\tpi1\thetVar1\tss1\tprivate1\tthetaH1\ttajd1\t\tH1\tHapCount1\tZnS1\tpi2h\tetVar2ss2\tprivate2\tthetaH2\ttajd2\tH2\tHapCount2\tZnS2\tFst\tsnn\tdxy_mean\tdxy_min\tgmin\tzx\tdd1\tdd2\tddRank1\tddRank2\tibsMaxB\tibsMean1\tibsMean2\tClass")

# Update stats-file with column on introgression-class
for line in Windows_file:
	data = line.split()
	window = data[0] + ":" + data[1] + "-" + data[2]

	if window in Introgression1: # Class 1
		output.write(line.rstrip())
		output.write("\t1\n")
	elif window in Introgression2: # Class 2
		output.write(line.rstrip())
		output.write("\t2\n")
	else:
		if "chrom" not in line: # Class 0
			output.write(line.rstrip())
			output.write("\t0\n")
 
Windows_file.close()
output.close()
```
Before I can plot the summary statistics in R, I need to change a few things.
```
# Replace . with , (R sees strings with a . as characters not numbers) - First sed-command
# Replace "inf" and "-nan" with NA - Second and third sed-commands
sed '-e s/\./,/g' -e 's/inf/NA/g' -e 's/-nan/NA/g' Window_Stats_Classification.txt > Temp.txt
rm Window_Stats_Classification.txt
mv Temp.txt Window_Stats_Classification.txt
```
Now, I can plot different summary stats in R to compare the introgression-classes. During the data-processing, I saved the Window_Stats_Classification-file in Excel-format. This is not necessary, it should also work with a standard txt-file.
```R
# Load libraries
library(ggplot2)
library(openxlsx)

# Load Data
All_Windows <- read.xlsx("Window_Stats_Classification.xlsx")
lapply(All_Windows, class) # Check if all the data are correct (should be numeric)
All_Windows$Class <- as.factor(All_Windows$Class)

# Remove introgression-class NA (=windows that could not be classified)
All_Windows <- All_Windows[ which(All_Windows$Class != 'NA'),]

# Plot summary statistic of your choice
ggplot(All_Windows, aes(x=Class, y=dxy_mean)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Introgression Class")
```
If necessary, I can also test for significant differences between the introgression classes with a Kruskal-Wallis and a post-hoc Tukey test.
```R
# Library for Tukey test
library(DescTools)

# Statistical tests
kruskal.test(All_Windows$dxy_mean ~ All_Windows$Class)
NemenyiTest(x = All_Windows$dxy_mean, g = All_Windows$Class, dist="tukey")
```

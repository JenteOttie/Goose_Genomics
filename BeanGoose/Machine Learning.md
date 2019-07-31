# Finding introgression regions using a Machine Learning approach
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
1.Introgression from population 1 to 2 
2.Introgression from population 2 to 1 
3.No Introgression
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

### Checking the results
To see which summary statistics were most important in classifying the windows, I run the python-script [getFeatureRankings.py](https://github.com/kern-lab/FILET/blob/master/getFeatureRankings.py) (written by Dan Schrider) on the file getFeatureRankings.py. This outputs a table that ranks the summary stats and shows their contribution to the classification.

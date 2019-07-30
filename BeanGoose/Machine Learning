# Finding introgression regions using a Machine Learning approach
## Background
In the previous step - demographic analysis - I showed that Taiga and Tundra Bean Goose have exchanged DNA during their evolutionary history.
To explore which genomic regions have been exchanged, I will apply a machine learning approach: FILET (Finding Introgressed Loci via Extra-Trees).
This software combines information from a suite of population genetic summary statistics that capture the genetic variation between two species (Schrider et al., 2018).
The machine learning algorithm relies on an Extra-Trees classifier (Geurts, Ernst, & Wehenkel, 2006), an extension of the Random Forest approach (Breiman, 2001). 
Random Forest is a collection of learning techniques that creates a forest of randomly generated decision trees that attempt to classify data points into predetermined groups. 
In this case, these groups correspond to introgression from population 1 into population 2, introgression from population 2 into population 1, and no introgression.

Here is a link to the FILET-paper: https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1007341

## Simulating training datasets
The algorithm needs to be trained. I will use simulated data as a training set. The data is simulated using [msmove](https://github.com/geneva/msmove).
Three scenarios are simulated: 1.Introgression from population 1 to 2 2.Introgression from population 2 to 1 3.No Introgression
```
msmove 18 10000 -t 50 -I 2 9 9 -n 1 1.40176528 -n 2 0.6912292 -ma x 5.80351132 26.15397712 x -ev 0.14821583 1 2 1 >mig12.txt
msmove 18 10000 -t 50 -I 2 9 9 -n 1 1.40176528 -n 2 0.6912292 -ma x 5.80351132 26.15397712 x -ev 0.14821583 2 1 1 >mig21.txt
msmove 18 10000 -t 0.1 -I 2 9 9 0.0001 -n 1 1.40176528 -n 2 0.6912292 -ev 0.14821583 1 2 0 >Nomig.txt
```

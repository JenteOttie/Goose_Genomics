# Demographic Analyses
To infer the evolutionary history of Taiga and Tundra Bean Goose, I will do demographic modelling using [dadi](https://bitbucket.org/gutenkunstlab/dadi/src/master/).
This software models allele frequency changes using diffusion approximation. For more information about this approach, you can read this [paper](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1000695).

![Dadi](https://journals.plos.org/plosgenetics/article/figure/image?size=large&id=10.1371/journal.pgen.1000695.g001)

## Preparing data
Because demographic modelling assumes neutral evolution, it is advised to use non-coding SNPs. I will select these SNPs using the software [snpEff](http://snpeff.sourceforge.net/).

**Installing snpEff**
```
wget http://sourceforge.net/projects/snpeff/files/snpEff_latest_core.zip
unzip snpEff_latest_core.zip
```

**Create database for Swan Goose genome**

The software snpEff has several databases installed. You can check them with the command.
```
$ java -jar snpEff.jar databases
```
Because the Swan Goose (Anser cygnoides) genome is not installed, I need to create it myself. First, I adjust the file snpEff.config. I add the following lines in the database section:
```
# Swan Goose genome, version 1
ansCyg1.genome : Swan Goose
```
Then I make folders data/ and genome/, which will hold the data from the Swan Goose genome. In the data-folder I save the gff-file and in the genome-folder the fasta-file. I change the names of these files to genes.gff.gz and ansCyg1.fa.gz.
Download data from https://www.ncbi.nlm.nih.gov/genome/?term=txid8845[orgn]
```
$ mv GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.gff.gz genes.gff.gz
$ mv GCF_000971095.1_AnsCyg_PRJNA183603_v1.0_genomic.fna.gz ansCyg1.fa.gz
```
Finally, I can build the database for this genome
```
$ java -jar snpEff.jar build -gff3 -v ansCyg1
```

**Annotating SNPs**

Using the just created database, I can annotate the SNPs in the vcf-file.
```
java -Xmx4g -jar snpEff/snpEff.jar -v ansCyg1 mydata.vcf >test.annotated.vcf
```

**Filtering SNPs**

For the demographic modelling, I will only use non-coding SNPs. So, I will have to remove all SNPs in coding regions. This can be done with the grep comment for multiple strings. This command will remove all the coding SNPs.
```
grep -v '3_prime_UTR_variant\|5_prime_UTR_premature_start_codon_gain_variant\|5_prime_UTR_variant\|initiator_codon_variant\|intragenic_variant\|missense_variant\|non_coding_transcript_exon_variant\|splice_acceptor_variant\|splice_donor_variant\|splice_region_variant\|start_lost\|stop_gained\|stop_lost\|stop_retained_variant\|synonymous_variant' BeanGoose_Filtered_Annotated.vcf >BeanGoose_Filtered_NonCoding.vcf
```

The remaining 12,114,494 SNPs will be filtered based on Hardy-Weinberg equilibrium (hwe), minor allele frequency (maf) and linkage disequilibrium (LD) using vcftools and plink. Make sure you use the right input-files for each step (for simplicity I just used mydata and result in the code below).

```
# Only use biallelic SNPs
vcftools --vcf mydata.vcf --max-missing-count 0 --max-alleles 2 --min-alleles 2 --remove-filtered-all --plink --out result
# Convert to BED-format
plink --file mydata --make-bed --out mydata
# Filter on hwe (0.01) and maf (0.05)
plink --bfile mydata --hwe 0.01 --maf 0.05 --makebed --out mydata_hwe_maf
# Change first column
awk '$1="1"' mydata.bim > mydata.bim
# LD pruning
plink --bfile mydata --indep-pairwise 50 10 0.05 --make-bed --out result
# Keep passed SNPs
plink --bfile mydata --extract BeanGoose_Filtered_hwe_maf_LD.prune.in --make-bed --out result 
# Convert bed-files back to vcf
plink --bfile mydata --recode vcf --out result
# Convert vcf to input for dadi
perl Convert_vcf_to_dadi.pl mydata.vcf Populations.txt
```
This leads to a final dataset of 5,397,934 SNPs.

&nbsp;

## Demographic Modelling
Dadi is based on Python-scripts. For these analyses, I need to start a virtual environment.
```
source dadi/bin/activate
```
Here is an example of a Python-script to run a model of isolation with asymmetrical migration. I will test a range of different scenarios for the Bean Geese.
```python
# Python-script to test demographic models using dadi.
# Run this script in virtual environment

# Import necessary packages
import numpy
import dadi

# Load the dataset
# File created with Perl-script Convert_vcf_to_dadi.pl
dd = dadi.Misc.make_data_dict('/proj/sllstore2017033/nobackup/work/jente/Bean_Goose_Project/Modelling/Dadi/BeanGoose_Filtered_NonCoding_TwoAlleles_hwe_maf_LD_0.5.vcf.data')

# Convert dictionary into Spectrum.
# Projections indicates sample size per population
# Polarized=False indicates that SNPs are not polarized with an outgroup. (i.e. spectrum is folded)
fs = dadi.Spectrum.from_data_dict(dd, pop_ids=['Fabalis','Rossicus'],
                            projections=[9,9],
                            polarized=False)

print('Data successfully loaded and converted to spectrum')

ns = [9,9]

###########################################################
# MODEL B - Isolation with asymmetrical migration
print('Running Model B: Isolation with asymmetrical migration')

def model_B(params, ns, pts):
    nu1, nu2, Ts, mBG, mGB = params

    xx = dadi.Numerics.default_grid(pts)

    phi = dadi.PhiManip.phi_1D(xx)
    phi = dadi.PhiManip.phi_1D_to_2D(xx, phi)

    phi = dadi.Integration.two_pops(phi, xx, Ts, nu1, nu2, m12=mBG, m21=mGB)

    fs = dadi.Spectrum.from_phi(phi, ns, (xx,xx))
    return fs

# Set bounds to parameters: nu1, nu2, Ts, m12, m21
upper_bound = [100, 100, 50, 50, 50]
lower_bound = [0.1, 0.1, 0, 0, 0]
p0 = [1, 1, 10, 15, 15]

# Set grid size
pts_l = [40,50,60]

# Prepare function
func = model_B
func_ex = dadi.Numerics.make_extrap_log_func(func)

# Perturb our parameters before optimization.
p0 = dadi.Misc.perturb_params(p0, fold=1, upper_bound=upper_bound,
                              lower_bound=lower_bound)
# Run optimization process
print('Beginning optimization ************************************************')
popt = dadi.Inference.optimize_log(p0,fs , func_ex, pts_l,
                                   lower_bound=lower_bound,
                                   upper_bound=upper_bound,
                                   verbose=len(p0), maxiter=3)
print('Finshed optimization **************************************************')

print('Best-fit parameters: {0}'.format(popt))

model = func_ex(popt, ns, pts_l)
theta = dadi.Inference.optimal_sfs_scaling(model, fs)

ll_model_B = dadi.Inference.ll_multinom(model, fs)

print('LL Model B: '), round(ll_model_B,0)

############################################################################
```

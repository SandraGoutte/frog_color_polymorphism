# Frog color polymorphism

This repository contains scripts to conduct analyses on the evolution and genomics of the green/brown color polymorphism in frogs, as described in Goutte and Boissinot (in review). Below is a brief description of these scripts and how to use them.

## Content
* [Dependencies](#dependencies)
* [Data](#data)
* [Evolution of green coloration in anurans](#evolution-of-green-coloration-in-anurans)
* [Diversification analysis](#diversification-analysis)
* [Genome-wide Association study of dorsal coloration in _Ptychadena robeensis_](#genome-wide-association-study-of-dorsal-coloration-in-_Ptychadena-robeensis_)
* [FoxD3 expression analysis](#foxD3-expression-analysis)



## Dependencies
The scripts in this directory use the following software and assume they are installed:

* [betascan](https://github.com/ksiewert/BetaScan)
* [bcftools](http://www.htslib.org/)
* [bcftools](http://www.htslib.org/)
* [beagle](https://faculty.washington.edu/browning/beagle/beagle.html)
* [bgzip](http://www.htslib.org/)
* [glactools](https://github.com/grenaud/glactools)
* [pixy](https://pixy.readthedocs.io/en/latest/)
* [plink](https://www.cog-genomics.org/plink/)
* [R](https://cran.r-project.org/)
* [tabix](http://www.htslib.org/)
* [vcf-kit]
* [vcftools](https://vcftools.github.io/index.html)

## Data

* Anurans' molecular phylogeny from [Portik et al. 2023](https://www.sciencedirect.com/science/article/pii/S1055790323002075)
* `TableS1.csv` - contains the color and habitat data for anurans
* [Robe's grass frog (_Ptychadena robeensis_) reference genome assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_036250615.1/)
* `robeensis.vcf.gz` - contains low-coverage resequencing data for 100 _P. robeensis_ individuals
* `phenotype_robeensis.txt` - contains color data (green or brown) for the 100 _P. robeensis_ sequenced individuals
* `countData.csv` - contains raw count data from RNAseq experiment
* `phenotypes2.csv` - contains information on the color, developmental stage, vertebral stripe width and zone (inside or outside the stripe) for each sample included in the RNAseq experiment
* `highcov10_genotype.vcf` - contains higher-coverage whole-genome sequencing data for 10 
_P. robeensis_ individuals
* `erlangeri_nana_levenorum_robeensis_highcov_genotype_all_sites.vcf.gz`  - contains higher-coverage whole-genome sequencing data for 4 polymorphic _Ptychadena_ species


## Evolution of green coloration in anurans
This analysis uses the molecular phylogeny for anurans from [Portik et al. 2023](https://www.sciencedirect.com/science/article/pii/S1055790323002075) and the dataset `TableS1.csv`. 

1. `Comparative_analysis_green.R` fits and compares Mk models of evolution of green coloration in anurans, reconstructs ancestral states, and plots the number of changes between color states and the evolutionary time spent in each state.
2. `Comparative_analysis_habitat_green.R` fits and compare Mk models of evolution for habitat preferences in anurans, reconstructs the ancestral states of habitat preferences, and fits and compares models of habitat-green color joint evolution.

## Diversification analysis
This analysis uses the same datasets as the [Evolution of green coloration in anurans](#evolution-of-green-coloration-in-anurans). 

`Diversification_analysis.R` fits and compare color state-dependent and independent models of diversification.

## Genome-wide Association study of dorsal coloration in _Ptychadena robeensis_

1. `GWAS.sh` filters the VCF file `robeensis.vcf.gz` and runs a genome-wide association study on _Ptychadena robeensis_ 100 individuals using the `phenotype_robeensis.txt` phenotype file.
2. `Plot_GWAS.R` imports the output of the GWAS and generates a Manhattan plot. 

## FoxD3 expression analysis

`Foxd3_expression_analysis` analyses the raw counts data from the RNAseq experiment and `phenotype2.csv` dataset. The  script normalizes the data and compares FoxD3 expression levels between green and brown frog skin across developmental stages, as well as between outside and inside the vertebral stripe in wide-striped individuals. 

## Selection analysis

This analysis is run on `highcov10_genotype.vcf` and `erlangeri_nana_levenorum_robeensis_highcov_genotype_all_sites.vcf.gz`

1. `Sliding_windows_stats` calculates Pi, Dxy, and Fst for 4 _Ptychadena_ species that share the same green/brown color polymorphism, and Tajima's D for _P. robeensis_ on sliding windows along the genome
2. `Run_Betascan.sh` runs Betascan on the _P. robeensis_
3. `Plot_stats_sliding_windows` loads all statistics computed on sliding windows and creates a mutli-panel plot for _P. robeensis_ and for the 3 other polymorphic _Ptychadena_ species

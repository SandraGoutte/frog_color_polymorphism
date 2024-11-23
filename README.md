# Frog color polymorphism

This repository contains scripts to conduct analyses on the evolution and genomics of the green/brown color polymorphism in frogs, as described in Goutte and Boissinot (in review). Below is a brief description of these scripts and how to use them. Shell and R scripts are in the [shell](https://github.com/SandraGoutte/frog_color_polymorphism/tree/main/shell) and [R](https://github.com/SandraGoutte/frog_color_polymorphism/tree/main/R) directories, respectively.

## Content
* [Dependencies](#dependencies)
* [Data used in this study](#data-used-in-this-study)
* [Evolution of the green coloration in anurans](#evolution-of-the-green-coloration-in-anurans)
* [Diversification analysis](#diversification-analysis)
* [Genome-wide Association Study of dorsal coloration in _Ptychadena robeensis_](#genome-wide-association-study-of-dorsal-coloration-in-_Ptychadena-robeensis_)
* [FoxD3 expression analysis](#foxD3-expression-analysis)
* [Selection analysis](#selection-analysis)
* [Recombination rate mapping](#recombination-rate-mapping)
* [TWISST analysis](#twisst-analysis)
* [Haplotype tree reconstruction](#haplotype-tree-reconstruction)



## Dependencies
The scripts in this directory use the following software and assume they are installed:

* [betascan](https://github.com/ksiewert/BetaScan)
* [bcftools](http://www.htslib.org/)
* [bcftools](http://www.htslib.org/)
* [beagle](https://faculty.washington.edu/browning/beagle/beagle.html)
* [bgzip](http://www.htslib.org/)
* [glactools](https://github.com/grenaud/glactools)
* [LDhat](https://github.com/auton1/LDhat)
* [pixy](https://pixy.readthedocs.io/en/latest/)
* [plink](https://www.cog-genomics.org/plink/)
* [R](https://cran.r-project.org/)
* [raxml-ng](https://github.com/amkozlov/raxml-ng)
* [tabix](http://www.htslib.org/)
* [twisst](https://github.com/simonhmartin/twisst)
* [vcf2fasta](https://github.com/santiagosnchez/vcf2fasta)
* [vcf-kit](https://github.com/AndersenLab/VCF-kit)
* [vcftools](https://vcftools.github.io/index.html)

## Data used in this study

* Anurans' molecular phylogeny from [Portik et al. 2023](https://www.sciencedirect.com/science/article/pii/S1055790323002075)
* `TableS1.csv` - contains the color and habitat data for anurans
* [Robe's grass frog (_Ptychadena robeensis_) reference genome assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_036250615.1/)
* `robeensis.vcf.gz` - contains low-coverage resequencing data for 100 _P. robeensis_ individuals
* `phenotype_robeensis.txt` - contains color data (green/brown) for the 100 _P. robeensis_ sequenced individuals
* `countData.csv` - contains raw count data from RNAseq experiment
* `phenotypes2.csv` - contains information on the color, developmental stage, vertebral stripe width and zone (inside or outside the stripe) for each sample included in the RNAseq experiment
* `highcov10_genotype.vcf` - contains higher-coverage whole-genome sequencing data for 10 
_P. robeensis_ individuals
* `4sp_highcov_all_sites.vcf.gz` - contains higher-coverage whole-genome sequencing data for 4 polymorphic _Ptychadena_ species, including invariant sites (necessary for Pixy)
* `highcov_12sp.vcf.gz` - contains higher-coverage whole-genome sequencing data for all 12 species of the Ethiopian Highlands _Ptychadena_ radiation
* `ptychadena_4sp.vcf` - contains higher-coverage whole-genome sequencing data for 4 polymorphic _Ptychadena_ species (no invariant sites)


## Evolution of the green coloration in anurans
This analysis uses the molecular phylogeny for anurans from [Portik et al. 2023](https://www.sciencedirect.com/science/article/pii/S1055790323002075) and the dataset `TableS1.csv`. 

1. `Comparative_analysis_green.R` fits and compares Mk models of evolution of the green coloration in anurans, reconstructs ancestral states, and plots the number of changes between color states and the evolutionary time spent in each state.
2. `Comparative_analysis_habitat_green.R` fits and compare Mk models of evolution for habitat preferences in anurans, reconstructs the ancestral states of habitat preferences, and fits and compares models of habitat-green color joint evolution.

## Diversification analysis
This analysis uses the molecular phylogeny for anurans from [Portik et al. 2023](https://www.sciencedirect.com/science/article/pii/S1055790323002075) and the dataset `TableS1.csv`. 

`Diversification_analysis.R` fits and compare color state-dependent and independent models of diversification.

## Genome-wide Association Study of dorsal coloration in _Ptychadena robeensis_

This analysis uses the `robeensis.vcf.gz` and `phenotype_robeensis.txt` files.

1. `GWAS.sh` filters the VCF file and runs a genome-wide association study on the dorsal coloration of 100 _Ptychadena robeensis_ individuals.
2. `Plot_GWAS.R` imports the output of the GWAS and generates a Manhattan plot. 

## FoxD3 expression analysis

This analysis uses the `countData.csv` and `phenotype2.csv` datasets.

`Foxd3_expression_analysis` normalizes the data and compares FoxD3 expression levels between green and brown frog skin across developmental stages, as well as between outside and inside the vertebral stripe in wide-striped individuals. 

## Selection analysis

This analysis is run on `4sp_highcov_all_sites.vcf.gz` and `highcov10_genotype.vcf`.

1. `Create_bed_file_sliding_window.R` creates a BED file with coordinates for a custom 3kb sliding window with 1kb overlap between adjacent windows to use in Pixy.
2. `Sliding_windows_stats.sh` calculates Pi, Dxy, and Fst for 4 _Ptychadena_ species that share the same green/brown color polymorphism, and Tajima's D for _P. robeensis_ on sliding windows along the genome.
3. `robeensis_thetamap.R` Calculates a mutation map for _P. robeensis_ to be used in Betascan.
4. `Run_Betascan.sh` runs Betascan on the _P. robeensis_ genome.
5. `plot_Betascores.R` imports the ouput from Betascan, calculates the first percentile of betascore values across the genome, and plots betascores.
6. `Plot_stats_sliding_windows.R` loads all statistics computed on sliding windows and creates a mutli-panel plot for _P. robeensis_ and for the 3 other polymorphic _Ptychadena_ species.

## Recombination rate mapping

This analysis is run on `highcov10_genotype.vcf`.

1. `Recombination_rate.sh` runs LDhat on the _P. robeensis_ genome and outputs a recombination map file per chromosome. This scipt uses `get_vcf_per_chr.sh`, `split_vcf_all_chr.sh`, `run_interval_all_chr.sh`, `joint_rho_SNPposition_all_chr.sh`, `remove_2_first_lines.sh`, `trim_result.sh`, `cat_all_results_one_file.sh`.
2. `plot_recombination_map.R` plots the recombination map in the region of interest.

## TWISST analysis

This analysis uses `ptychadena_4sp.vcf`.

1. `Run_TWISST_on_sliding_windows.sh` phases the input VCF and runs TWISST on the 4 polymrophic _Ptychadena_ species (including green and brown individuals for each species). In order to make the windows artificially overlap (with 1 SNP step), the script creates 50 versions of the VCF starting at 1 SNP intervall from each other and TWISST is run on 50 SNPs non-overlapping windows on each VCF, results are combined at the plotting step.
2. `Plot_TWISST_output.R` imports the output from TWISST, identifies tree topologies in which green individuals of all 4 species group together, sums the weights of green-grouping topologies and plot the result.

## Haplotype tree reconstruction

This analysis uses `highcov_12sp.vcf.gz`.

`Reconstruct_haplotype_tree.sh` phases the VCF including 12 Ethiopian Highlands _Ptychadena_ species, selects the region of interest and splits the haplotypes before reconstructing a phylogenetic tree of the haplotypes.

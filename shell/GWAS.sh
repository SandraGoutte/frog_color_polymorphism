#####################################################################
#####
#####               GWAS color Ptychadena robeensis 2024
#####
#####################################################################
### What this script does: 
### 1. filters the VCF 
### 2. convert the VCF into plink format
### 3. run a GWAS on color trait using plink

### Dependencies:
### vcftools
### plink

### Data
## robeensis.vcf.gz contains low-coverage resequencing data for 100 P. robeensis individuals
## check the average coverage depth per sample and site quality in our VCF file
vcftools --gzvcf robeensis.vcf.gz --depth --site-quality --out robeensis

## check the number of sites in the VCF file
zcat robeensis.vcf.gz | grep -v '#' | wc -l
#48608304

## phenotype_robeensis.txt has the phenotype data associated with our VCF file
## it contains 1 row/individual and 1 column/trait (color: 1 = brown, 2 = green)

## check our phenotype data
head phenotype_robeensis.txt
#FID	#IID	color	
11_933_S57	11_933_S57	1	
11_934_S77	11_934_S77	2	
13_106_S191	13_106_S191	2	

#####################################################################
### 1. Filter VCF file
#####################################################################

## filter out non biallelic sites, sites with a quality <30, sites with > 40% missing data,
## and sites with a coverage > 15x.
vcftools --gzvcf robeensis.vcf.gz \
        --max-alleles 2 \
        --max-meanDP 15 \
        --min-alleles 2 \
        --max-missing 0.6 \
        --minQ 30 \
        --out robeensis_filtered \
        --recode

## check the number of sites in the VCF file
cat robeensis_recode.vcf | grep -v '#' | wc -l
#6700191

## compress the VCF file
bgzip -c  robeensis_filtered.recode.vcf > robeensis_filtered_recode.vcf.gz 

#####################################################################
### 2. Convert vcf to plink format
#####################################################################
./plink --vcf robeensis_filtered_recode.vcf.gz --make-bed --out  robeensis_filtered --allow-extra-chr --double-id --allow-no-sex

## check for relatedness in our data 
./plink --bfile robeensis_filtered --make-rel --allow-extra-chr --out robeensis_filtered_rel
### all individuals are lower than 0.2, which is our threshold

## check the distribution of HWE p-values of all SNPs.
./plink --bfile robeensis_filtered --hardy --allow-extra-chr

## filters out only SNPs which deviate extremely from HWE. 
./plink --bfile robeensis_filtered --hwe 1e-20 --make-bed --out robeensis_filtered_hwe --allow-extra-chr 
#3708 variants removed due to Hardy-Weinberg exact test 

#####################################################################
### 3. Run GWAS on the filtered data
#####################################################################
./plink --bfile robeensis_filtered_hwe --assoc --pheno phenotype_robeensis.txt --adjust --all-pheno --out robeensis_filtered_hwe --allow-extra-chr --allow-no-sex 

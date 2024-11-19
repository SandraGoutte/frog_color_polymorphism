#####################################################################
##### 
#####     Calculating Pi, Dxy, Fst, Tajima's D on sliding windows
##### 
#####################################################################
### What this script does: 
### 1. Install Pixy
### 2. Filter VCF and run Pixy
### 3. Get Tajima'D in sliding windows using VCF-KIT

### Dependencies:
### vcftools
### Pixy
### vcf-kit

### Data:
### erlangeri_nana_levenorum_robeensis_highcov_genotype_all_sites.vcf.gz 
### is a VCF file including our 4 polymoprhic Ptychadena species 
### and containing invariant sites (necessary for Pixy)

#####################################################################
### 1. Install Pixy 
#####################################################################
## create a dedicated environment
conda create --name pixy
conda activate pixy

## install Pixy
conda install --yes -c conda-forge pixy
conda install --yes -c bioconda htslib

#####################################################################
### 2. Run Pixy on VCF
#####################################################################
###  VCF NEEDS TO INCLUDE NON-VARIANT SITES  ###

#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --mem 500GB

## load the environment within the bash job
source  activate /.../envs/pixy

module load gencore
module load gencore_variant_detection/1.0

## filter our VCF: remove sites with > 60% missing, with depth > 15 , or with more that 2 alleles
vcftools --gzvcf erlangeri_nana_levenorum_robeensis_highcov_genotype_all_sites.vcf.gz --recode --recode-INFO-all --mac 3 --max-alleles 2 --max-missing 0.4 --max-meanDP 15 --out highcov_4sp_allsites_masmissing04_maxDP15

## extract only chromosome 7
vcftools --gzvcf erlangeri_nana_levenorum_robeensis_highcov_genotype_all_sites.vcf.gz --chr Scaffold_4__1_contigs__length_105873206 --recode --recode-INFO-all --mac 3 --max-alleles 2 --max-missing 0.4 --max-meanDP 15 --out highcov_4sp_chr7_allsites_masmissing04_maxDP15

## bgzip the VCFs
bgzip  highcov_4sp_allsites_masmissing04_maxDP15.recode.vcf
tabix  highcov_4sp_allsites_masmissing04_maxDP15.recode.vcf.gz

bgzip  highcov_4sp_chr7_allsites_masmissing04_maxDP15.recode.vcf
tabix  highcov_4sp_chr7_allsites_masmissing04_maxDP15.recode.vcf.gz

## run Pixy using a custom-made sliding window of 3kb overlapping 1kb
## species and color on chromosome 7
pixy --stats pi fst dxy \
--vcf highcov_4sp_chr7_allsites_masmissing04_maxDP15.recode.vcf.gz \
--populations pops_4sp_color.txt \
--bed_file sliding_windows_roi_step1000.bed \
--n_cores 4 \
--bypass_invariant_check 'yes' \
--chromosomes "Scaffold_4__1_contigs__length_105873206" \
--output_prefix pixy4spallsitescustomwindowchr7colorstep1000_filtered

## species (both colors together) on chromosome 7
pixy --stats pi fst dxy \
--vcf highcov_4sp_chr7_allsites_masmissing04_maxDP15.recode.vcf.gz \
--populations pops_4sp.txt \
--bed_file sliding_windows_roi_step1000.bed \
--n_cores 4 \
--bypass_invariant_check 'yes' \
--chromosomes "Scaffold_4__1_contigs__length_105873206" \
--output_prefix pixy4spallsitescustomwindowchr7colorstep1000_filtered_species

## run Pixy on the entire genome to get average values 
pixy --stats pi fst dxy \
--vcf highcov_4sp_chr7_allsites_masmissing04_maxDP15.recode.vcf.gz \
--populations pops_4sp.txt \
--bed_file sliding_windows_fullgenome_w3000_step1000.bed \
--n_cores 4 \
--bypass_invariant_check 'yes' \
--output_prefix pixy4spallsitescustomwindow_w3000_step1000_full_genome_species

#####################################################################
### 3. Get Tajima'D in sliding windows using VCF-KIT
#####################################################################
## install VCF-KIT 
pip install numpy # You may need to install numpy independently.
pip install VCF-kit
vk setup

## get Tajima's D on a 3kb sliding window (1kb overlap) on chromosome 7
vk tajima 3,000 1,000 highcov_robeensis_chr7.recode.vcf > TajimasD_robeensis_highcov_chr7_2.txt

## get Tajima's D on a 3kb sliding window (1kb overlap) over the entire genome
vk tajima 3,000 1,000 highcov10_genotype.vcf > TajimasD_robeensis_highcov.txt

## To get the Tajima's D for each color morph, lets subsample our VCF into green and brown ones
## green
vcftools --gzvcf highcov_4sp_chr7_allsites_masmissing04_maxDP15.recode.vcf.gz --chr Scaffold_4__1_contigs__length_105873206 --recode --recode-INFO-all --min-alleles 2 --max-alleles 2 --max-missing 0.4 --max-meanDP 15 --keep green_robeensis_highcov.txt --out highcov_robeensis_chr7_green

## brown
vcftools --gzvcf highcov_4sp_chr7_allsites_masmissing04_maxDP15.recode.vcf.gz --chr Scaffold_4__1_contigs__length_105873206 --recode --recode-INFO-all --min-alleles 2 --max-alleles 2 --max-missing 0.4 --max-meanDP 15 --keep brown_robeensis_highcov.txt --out highcov_robeensis_chr7_brown

## get Tajima's D on a 3kb sliding window (1kb overlap) on chromosome 7 (brown)
vk tajima 3,000 1,000 highcov_robeensis_chr7_brown.recode.vcf > TajimasD_robeensis_highcov_chr7_brown.txt

## get Tajima's D on a 3kb sliding window (1kb overlap) on chromosome 7 (green)
vk tajima 3,000 1,000 highcov_robeensis_chr7_green.recode.vcf > TajimasD_robeensis_highcov_chr7_green.txt


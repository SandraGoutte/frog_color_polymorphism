#####################################################################
#####  
#####       Recombination rate analysis for Ptychadena robeensis 
##### 
#####################################################################
### What this script does: 
### 1. Install LDhat
### 2. Prepare VCF files to run LDhat
### 3. Prepare lookup table to run LDhat
### 4. Run LDhat
### 5. Joins all results into a single output per chromosome

### Dependencies:
### vcftools
### gcc
### LDhat (installation at step 1)

### Data
### highcov10_genotype.vcf : VCF file with 10 high coverage robeensis

#####################################################################
### 1. Install LDhat
#####################################################################
git clone https://github.com/auton1/LDhat
cd LDhat
make

#####################################################################
### 2. VCF file preparation for LDhat
#####################################################################
## check the average coverage to know what filter to use
vcftools --vcf highcov10_genotype.vcf --site-mean-depth
# check the output out.ldepth.out file in R (LDhat_plot.R)
# quartiles of the mean depth per site= 2.0  6.6  7.9  9.7 14.3
# use the filters: min depth = 6, max depth = 14

## filter VCF file: min depth = 6, max depth=14, max missing data in 50% individuals, keep only biallelic sites
vcftools --vcf highcov10_genotype.vcf --max-alleles 2 --min-alleles 2 --min-meanDP 6 --max-meanDP 14 --max-missing 0.5 --out highcov10_genotype_filtered --recode

## get a VCF file for each chr, each in its own directory
sbatch get_vcf_per_chr.sh

## split each chromosome in ~2000 SNPs-fragments 
# we use 450,000 bp fragments with 1000 bp overlap
# create the fragments by running the script split_vc_all_chr.sh
sbatch split_vcf_all_chr.sh

#####################################################################
### 3. Lookup table preparation for LDhat
#####################################################################
## get the lookup table for 10 individuals
## download the lookup table with the closest parameters from our dataset on the github here: https://github.com/auton1/LDhat/tree/master/lk_files
wget https://github.com/auton1/LDhat/blob/master/lk_files/lk_n50_t0.001.gz

## unzip the lookup table
gzip -d lk_n50_t0.001.gz

### generate a new lookup table with the premade one and our file
## Convert VCF file to LDhat format # requires a chromosome for this format, take chr1
vcftools --vcf /scratch/sg5533/NEW_GENOME_ANALYSES/ROBEENSIS_GENOME/HIGHCOV/highcov10_genotype_filtered.recode.vcf --out highcov10_filtered_chr1 --ldhat-geno --chr Scaffold_2__2_contigs__length_236349107

## use the convert function from LDhat to get a site.txt file
./LDhat/convert -seq highcov10_filtered_chr1.ldhat.sites -loc highcov10_filtered_chr1.ldhat.locs

## create a new lookup table from that premade one # for n=2x  nb of individuals because our species is diploid = 2x 10 = 20 sequences
./LDhat/lkgen -lk lk_n50_t0.001.txt -seq sites.txt

#####################################################################
### 4. Run LDhat on each fragment
#####################################################################
## run the LDhat interval function on all chromosomes
sbatch run_interval_all_chr.sh

## summarize the results from interval 
## 2 chromosomes at a time otherwise it takes too long
sbatch run_stat_interval_chr1-2.sh
sbatch run_stat_interval_chr3-4.sh
sbatch run_stat_interval_chr5-6.sh
sbatch run_stat_interval_chr7-8.sh
sbatch run_stat_interval_chr9-10.sh
sbatch run_stat_interval_chr11-12.sh

#####################################################################
### 5. Join all the results for plotting in R
#####################################################################
## the res.text file outputted by Stats is a value of rho per SNPs. We need to make a correspondence between the SNP position (given in the .locs file, in kb) and the rho in the res file.
sbatch joint_rho_SNPposition_all_chr.sh

## remove the 2 first line for each joint file (the header and an aberrant value of rho)
sbatch remove_2_first_lines.sh

## our files overlap 1000 bp in total, but that’s different numbers of lines for each file, depending on the density of SNPs. So we need to remove 500 bp on each side of each file (except the first and last) and multiply the position by 1000 to have bp instead of kb
sbatch trim_result.sh

## row bind all the tables to have a single file for each chr
sbatch cat_all_results_one_file.sh
# all end up in the first file chr11_highcov_1.joint.txt
# download and rename chr11_highcov_1-127.joint.txt and plot in R 

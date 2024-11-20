#####################################################################
#####
##### Test for balancing selection on the GCP locus in P. robeensis
#####
#####################################################################
### What this script does: 
### 1. Install Betascan
### 2. Filter and subsample obtain 1 VCF file per chromosome
### 3. Convert the VCF file into the BetaScan format
### 4. Get the mutation rate (theta) map for each chromosome
### 5. Run Betascan on the 12 chromososomes

### Dependencies:
### numpy
### vcftools
### glactools

### Data:
### highcov10_genotype.vcf contains "high-coverage" resequencing data 
### for 10 Ptychadena robeensis individuals

#####################################################################
### 1. Install Betascan
#####################################################################
## create a dedicated environement
conda create -n Betascan
conda activate Betascan

## install numpy if it isn't installed yet
conda install numpy

## download the BetaScan python file
wget https://github.com/ksiewert/BetaScan/blob/master/BetaScan.py

#####################################################################
### 2. Filter and subsample obtain 1 VCF file per chromosome
#####################################################################
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --mem 500GB

mkdir robeensis_all_chr

cd robeensis_all_chr
vcftools --vcf highcov10_genotype.vcf \
--chr Scaffold_2__2_contigs__length_236349107 \
--recode \
--mac 3 \
--max-missing 0.2 \
--max-meanDP 25 \
--out robeensis_chr1;

vcftools --vcf highcov10_genotype.vcf \
--chr Scaffold_1__1_contigs__length_194224858 \
--recode \
--mac 3 \
--max-missing 0.2 \
--max-meanDP 25 \
--out robeensis_chr2;

vcftools --vcf highcov10_genotype.vcf \
--chr Scaffold_5__2_contigs__length_177367226 \
--recode \
--mac 3 \
--max-missing 0.2 \
--max-meanDP 25 \
--out robeensis_chr3;

vcftools --vcf highcov10_genotype.vcf \
--chr Scaffold_3__2_contigs__length_170906183 \
--recode \
--mac 3 \
--max-missing 0.2 \
--max-meanDP 25 \
--out robeensis_chr4;

vcftools --vcf highcov10_genotype.vcf \
--chr Scaffold_6__2_contigs__length_157086139 \
--recode \
--mac 3 \
--max-missing 0.2 \
--max-meanDP 25 \
--out robeensis_chr5;

vcftools --vcf highcov10_genotype.vcf \
--chr Scaffold_7__2_contigs__length_128225229 \
--recode \
--mac 3 \
--max-missing 0.2 \
--max-meanDP 25 \
--out robeensis_chr6; 

vcftools --vcf highcov10_genotype.vcf \
--chr Scaffold_4__1_contigs__length_105873206 \
--recode \
--mac 3 \
--max-missing 0.2 \
--max-meanDP 25 \
--out robeensis_chr7;

vcftools --vcf highcov10_genotype.vcf \
--chr Scaffold_10__2_contigs__length_97064778 \
--recode \
--mac 3 \
--max-missing 0.2 \
--max-meanDP 25 \
--out robeensis_chr8;

vcftools --vcf highcov10_genotype.vcf \
--chr Scaffold_8__1_contigs__length_80766577 \
--recode \
--mac 3 \
--max-missing 0.2 \
--max-meanDP 25 \
--out robeensis_chr9;

vcftools --vcf highcov10_genotype.vcf \
--chr Scaffold_9__1_contigs__length_67344957 \
--recode \
--mac 3 \
--max-missing 0.2 \
--max-meanDP 25 \
--out robeensis_chr10;

vcftools --vcf highcov10_genotype.vcf \
--chr Scaffold_12__2_contigs__length_57092889 \
--recode \
--mac 3 \
--max-missing 0.2 \
--max-meanDP 25 \
--out robeensis_chr11;

vcftools --vcf highcov10_genotype.vcf \
--chr Scaffold_11__1_contigs__length_50805362 \
--recode \
--mac 3 \
--max-missing 0.2 \
--max-meanDP 25 \
--out robeensis_chr12;

#####################################################################
### 3. Convert the VCF file into the BetaScan format
#####################################################################
## install glactools
module load gcc/9.2.0 
git clone --depth 1 https://github.com/grenaud/glactools.git
cd glactools
make
cd ..

## format all the chromosomes VCFs for Betascan
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --mem 118GB

for i in {1..12}
do   
    ./glactools/glactools vcfm2acf --onlyGT --fai /.../ptychadena_dovetail_2.fa.fai /.../robeensis_all_chr/robeensis_chr$i.recode.vcf > /.../robeensis_all_chr/robeensis_chr$i.recode.vcf.acf.gz ;

    ./glactools/glactools acf2betascan --fold /.../robeensis_all_chr/robeensis_chr$i.recode.vcf.acf.gz | gzip > /.../robeensis_all_chr/robeensis_chr$i.recode.beta.txt.gz 
done


#####################################################################
### 4. Get the mutation rate (theta) map for each chromosome
#####################################################################
## first get the SNP density for windows of 12kb using VCFtools
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --mem 118GB

vcftools --vcf robeensis_chr1.recode.vcf \
--out robeensis_chr1_12kb_window --SNPdensity 12000;
vcftools --vcf robeensis_chr2.recode.vcf \
--out robeensis_chr2_12kb_window --SNPdensity 12000;
vcftools --vcf robeensis_chr3.recode.vcf \
--out robeensis_chr3_12kb_window --SNPdensity 12000;
vcftools --vcf robeensis_chr4.recode.vcf \
--out robeensis_chr4_12kb_window --SNPdensity 12000;
vcftools --vcf robeensis_chr5.recode.vcf \
--out robeensis_chr5_12kb_window --SNPdensity 12000;
vcftools --vcf robeensis_chr6.recode.vcf \
--out robeensis_chr6_12kb_window --SNPdensity 12000;
vcftools --vcf robeensis_chr7.recode.vcf \
--out robeensis_chr7_12kb_window --SNPdensity 12000;
vcftools --vcf robeensis_chr8.recode.vcf \
--out robeensis_chr8_12kb_window --SNPdensity 12000;
vcftools --vcf robeensis_chr9.recode.vcf \
--out robeensis_chr9_12kb_window --SNPdensity 12000;
vcftools --vcf robeensis_chr10.recode.vcf \
--out robeensis_chr10_12kb_window --SNPdensity 12000;
vcftools --vcf robeensis_chr11.recode.vcf \
--out robeensis_chr11_12kb_window --SNPdensity 12000;
vcftools --vcf robeensis_chr12.recode.vcf \
--out robeensis_chr12_12kb_window --SNPdensity 12000

## Download and modify the tables in R (robeensis_thetamap.R)
## we obtain a tab-delimited text file for each chromosome

## remove first column and first row of each file
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --mem 118GB

for i in {1..12}
do  

    gawk -i inplace 'NR>1'  /.../Betascan/robeensis_all_chr/theta_map_chr$i.12kb_window.txt; 
    gawk -i inplace '{$1=""}1' /.../Betascan/robeensis_all_chr/theta_map_chr$i.12kb_window.txt
 
done

#####################################################################
### 5. Run Betascan on the 12 chromososomes
#####################################################################
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --mem 118GB

for i in {1..12}
do   
 python  BetaScan.py -i /.../Betascan/robeensis_all_chr/robeensis_chr$i.recode.beta.txt.gz -w 3000 -fold -m 0.15 -theta_map /.../Betascan/robeensis_all_chr/theta_map_chr$i.12kb_window.txt -o /.../Betascan/robeensis_all_chr/robeensis_chr$i.w3000_thetamap.betascores.txt
done
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --mem 118GB

module load vcftools

mkdir chr1
cd chr1

vcftools --vcf highcov10_genotype_filtered.recode.vcf \
--out robeensis_highcov_chr1 \
--chr Scaffold_2__2_contigs__length_236349107 \
--recode

cd ..
mkdir chr2
cd chr2

vcftools --vcf highcov10_genotype_filtered.recode.vcf \
--out robeensis_highcov_chr2 \
--chr Scaffold_1__1_contigs__length_194224858 \
--recode

cd ..
mkdir chr3
cd chr3

vcftools --vcf highcov10_genotype_filtered.recode.vcf \
--out robeensis_highcov_chr3 \
--chr Scaffold_5__2_contigs__length_177367226 \
--recode

cd ..
mkdir chr4
cd chr4

vcftools --vcf highcov10_genotype_filtered.recode.vcf \
--out robeensis_highcov_chr4 \
--chr Scaffold_3__2_contigs__length_170906183 \
--recode

cd ..
mkdir chr5
cd chr5

vcftools --vcf highcov10_genotype_filtered.recode.vcf \
--out robeensis_highcov_chr5 \
--chr Scaffold_6__2_contigs__length_157086139 \
--recode

cd ..
mkdir chr6
cd chr6

vcftools --vcf highcov10_genotype_filtered.recode.vcf \
--out robeensis_highcov_chr6 \
--chr Scaffold_7__2_contigs__length_128225229 \
--recode

cd ..
mkdir chr7
cd chr7

vcftools --vcf highcov10_genotype_filtered.recode.vcf \
--out robeensis_highcov_chr7 \
--chr Scaffold_4__1_contigs__length_105873206 \
--recode

cd ..
mkdir chr8
cd chr8

vcftools --vcf highcov10_genotype_filtered.recode.vcf \
--out robeensis_highcov_chr8 \
--chr Scaffold_10__2_contigs__length_97064778 \
--recode

cd ..
mkdir chr9
cd chr9

vcftools --vcf highcov10_genotype_filtered.recode.vcf \
--out robeensis_highcov_chr9 \
--chr Scaffold_8__1_contigs__length_80766577 \
--recode

cd ..
mkdir chr10
cd chr10

vcftools --vcf highcov10_genotype_filtered.recode.vcf \
--out robeensis_highcov_chr10 \
--chr Scaffold_9__1_contigs__length_67344957 \
--recode

cd ..
mkdir chr11
cd chr11

vcftools --vcf highcov10_genotype_filtered.recode.vcf \
--out robeensis_highcov_chr11 \
--chr Scaffold_12__2_contigs__length_57092889 
\--recode

cd ..
mkdir chr12
cd chr12

vcftools --vcf highcov10_genotype_filtered.recode.vcf \
--out robeensis_highcov_chr12 \
--chr Scaffold_11__1_contigs__length_50805362 \
--recode

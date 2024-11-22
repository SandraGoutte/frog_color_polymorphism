#####################################################################
##### 
#####      Running TWISST on overlapping 50 SNPs sliding windows       
#####
#####################################################################
### What this script does: 
### 1. Phase VCF with Beagle
### 2. Rename the samples in VCF
### 3. Create a group file for TWISST
### 4. Extract VCFs
### 5. Install TWISST
### 6. Infer trees for TWISST
### 7. run TWISST 

### Dependencies:
### Beagle
### vcftools
### bcftools
### parseVCF.py from https://github.com/simonhmartin/genomics_general/
### phyml_sliding_windows.py from https://github.com/simonhmartin/genomics_general/

### Data:
### ptychadena_4sp.vcf is a VCF file including 
### "high-coverage" sequencing data of our 4 polymorphic species: 
### P. nana, P. erlangeri, P. levenorum, P. robeensis

#####################################################################
### 1. Phase VCF using Beagle
#####################################################################
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=500GB
#SBATCH --time=96:00:00

## load modules
module load gencore
module load gencore_variant_detection

## run Beagle
java -Xmx300g -jar beagle.22Jul22.46e.jar gt=ptychadena_4sp.vcf out=phased_ptychadenas4sp.vcf.gz impute=true nthreads=20 window=40 overlap=2 gp=false

#####################################################################
### 2. Rename samples in VCF
#####################################################################
## in order to avoid any issues with TWISST, rename the samples to 
## numbers
## check the sample names
vcf-query -l phased_ptychadenas4sp.vcf.gz 

## make a tab-delimited text file with old and new names: rename_samples_ptychadenas4sp.txt
1	13_106_S37
2	SB495_S36
3	SB653_S38
...

## replace the sample names in the VCF file
bcftools reheader -s rename_samples_ptychadenas4sp.txt phased_ptychadenas4sp.vcf.gz  > renamed_phased_ptychadenas4sp.vcf.gz   

## check the new sample names (should be only numbers)
vcf-query -l renamed_phased_ptychadenas4sp.vcf.gz  

#####################################################################
### 3. Create a group file as input for TWISST
#####################################################################
## use the list of samples in rename_samples_ptychadenas4sp.txt
## and create a tab-deliminted text file with 2 columns, no header
## for each sample (column 1), since our VCF is phased, create 2 haplotypes
## for example sample 1 becomes 1_A and 1_B
## for each haplotype, assign a group for TWISST (column 2) 
## in our case these correspond to species_color
head groups_TWISST_4sp.txt 
#1_A	erlangeri_brown
#1_B	erlangeri_brown
#10_A	nana_brown
#10_B	nana_brown

#####################################################################
#### 4. Extract 50 VCFs each starting at one of 50 first SNPs
#####################################################################
## To obtain 50 SNPs windows that overalp on 49 SNPs, we create 
## our own custom window system: I have a list of the SNPs position here: 
## phased_ptychadena_4sp_chr7_26000000_27000000.ldepth
## in R: 
## depth[1:50,]$POS
## write the data to a tab-delimited text file
## write.table(depth[1:50,]$POS, file = "50_first_snps.txt", 
## sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)

## extract VCFs starting at each of these SNPs (50 VCF files in total)
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --mem 400GB

## load modules
module load gencore
module load gencore_variant_detection/1.0

## create a new directory and place yourself in it
mkdir vcfs
cd vcfs

## path to the text file containing the positions of the 50 first SNPs
input_file="/.../50_first_snps.txt"

## path to our phased VCF file
vcf_file="renamed_phased_ptychadenas4sp.vcf.gz"
chromosome="Scaffold_4__1_contigs__length_105873206"

## read start and end positions from the input file into arrays
mapfile -t start_positions < <(awk '{print $1}' $input_file)

## loop through the arrays
for i in "${!start_positions[@]}"; do
    start="${start_positions[$i]}"
    end=27000000

    ## extract SNPs for the specified region
    vcftools --gzvcf $vcf_file \
    --recode-INFO-all \
    --recode \
    --chr $chromosome \
    --from-bp $start \
    --to-bp $end \
    --out renamed_phased_ptychadena_4sp_chr7_$i
done

## bgzip each vcf file
bgzip  renamed_phased_ptychadena_4sp_chr7_*.recode.vcf
tabix  renamed_phased_ptychadena_4sp_chr7_*.recode.vcf.gz

#####################################################################
### 5. Install TWISST
#####################################################################
## download the twisst package zip file from GitHub
wget https://github.com/simonhmartin/twisst/archive/v0.2.tar.gz

## extract the files from the zipped file
tar -xzf v0.2.tar.gz

## remove the zipped file
rm v0.2.tar.gz

#####################################################################
### 6. Infer trees to use in TWISST
#####################################################################
### Infer trees for each window
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --mem 400GB

for vcf in {0..49}; do
echo "Inferring trees with window size 50 for vcf $vcf"

python parseVCF.py -i /.../renamed_phased_ptychadena_4sp_chr7_$vcf.recode.vcf.gz | 
python phyml_sliding_windows.py \
--prefix renamed_phased_ptychadena_4sp_chr7_w50SNPs_$vcf.phyml_bionj \
--windType sites -w 50  --model GTR --optimise n
done

#####################################################################
### 7. run TWISST 
#####################################################################
## run TWISST on each VCF file with a 50 snp window
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --mem 400GB

for vcf in {0..49}; do
## run TWISST with 8 groups : P. robeensis, P. nana, P. levenorum, P. erlangeri. 
## for each species: green / brown
python twisst.py \
-t renamed_phased_ptychadena_4sp_chr7_w50SNPs_$vcf.phyml_bionj.trees.gz \
-w renamed_phased_ptychadena_4sp_chr7_w50SNPs_$vcf.phyml_bionj.weights.tsv \
-g robeensis_brown \
-g robeensis_green \
-g nana_brown \
-g nana_green \
-g levenorum_brown \
-g levenorum_green \
-g erlangeri_brown \
-g erlangeri_green \
--groupsFile groups_TWISST_4sp.txt
done

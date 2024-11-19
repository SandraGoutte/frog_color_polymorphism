#####################################################################
#####
#####   Reconstruct a haplotype tree for the region of interest     
#####
#####################################################################
### What this script does: 
### 1. Phase VCF with Beagle
### 2. Get a VCF file with only the region of interest
### 3. Split the haplotypes
### 4. convert the VCF file into fasta
### 5. Reconstruct a phylogeny for the region of interest

### Dependencies:
### Beagle
### vcftools
### vcf2fasta from https://github.com/santiagosnchez (install at step 4)
### raxml-ng (install at step 5)

### Data:
### highcov_12sp_genotype.recode.vcf.gz is a VCF including "high-coverage" 
### sequencing data from all 12 species of the Ethiopian Highlands Ptychadena
### radiation

#####################################################################
### 1. Phase VCF using Beagle
#####################################################################
java -Xmx300g -jar beagle.22Jul22.46e.jar gt=highcov_12sp_genotype.recode.vcf.gz out=phased_ptychadenas12sp.vcf.gz impute=true nthreads=20 window=40 overlap=2 gp=false

#####################################################################
### 2. Get a VCF file with only the region of interest
#####################################################################
vcftools --gvcf phased_ptychadenas12sp.vcf.gz --recode-INFO-all --recode --chr Scaffold_4__1_contigs__length_105873206 --from-bp 26414430 --to-bp 26417815 --out phased_12sp_chr7_win1_and_2_26414430-26417815
#kept 762 out of a possible 1582 Sites

#####################################################################
### 3. Split the haplotypes
#####################################################################
sed "s:|:\t:g" phased_12sp_chr7_win1_and_2_26414430-26417815_recode.vcf > phased_12sp_chr7_win1_and_2_26414430-26417815_recode.split.vcf

## check the resulting vcf file
less phased_12sp_chr7_win1_and_2_26414430-26417815_recode.split.vcf

## get the header of the file and then edit it to have the haplotypes A/B 
## edit manually
grep "^#CHROM" phased_12sp_chr7_win1_and_2_26414430-26417815_recode.split.vcf > header_12sp.txt

head header_12sp.txt
##CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	13-290_S4	13_106_S37	15-104	15-139_S1	15-158_S3	15-173_S2	15-260_S5	15-283_S4	15-328_S6	15-449_S5	16-126_S2	SB495_S36	

## edit to have each sample followed by _A and _B
#NOTE: Following symbols are not allowed in taxa names to ensure Newick compatibility:
#NOTE: " " (space), ";" (semicolon), ":" (colon), "," (comma), "()" (parentheses), "'" (quote)
## If there are spaces in between the names instead of tabs, edit in TextEdit the line with the names to remplace by tabs BEFORE RUNNING THE VCF2FASTA

head header_12sp.txt
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	13-290_S4_A 13-290_S4_B	13_106_S37_A	13_106_S37_B 15-104_A 15-104_B	15-139_S1_A	15-139_S1_B 15-158_S3_A	15-158_S3_B 15-173_S2_A

#####################################################################
### 4. convert the VCF file into fasta
#####################################################################
## compress and index VCF file
bgzip phased_12sp_chr7_win1_and_2_26414430-26417815_recode.split.vcf
tabix phased_12sp_chr7_win1_and_2_26414430-26417815_recode.split.vcf.gz

## create a dummy gff file as imput for vcf2fasta: dummy_gff.gff 
##gff-version 3
Scaffold_4__1_contigs__length_105873206	.	CDS	26414430	26417815	.	+	.	ID=mrna0001;Name=window12
Scaffold_4__1_contigs__length_105873206	.	CDS	26414430	26415921	.	+	.	ID=mrna0001;Name=window1
Scaffold_4__1_contigs__length_105873206	.	CDS	26416840	26417815	.	+	.	ID=mrna0001;Name=window2

## install the vcf2fasta python script
git clone https://github.com/santiagosnchez/vcf2fasta
cd vcf2fasta
pip3 install pysam art

## convert the vcf file into fasta
python vcf2fasta.py --fasta /.../ptychadena_dovetail_2.fa --vcf phased_12sp_chr7_win1_and_2_26414430-26417815_recode.split.vcf.gz --gff dummy_gff.gff --feat CDS

#####################################################################
### 5. Reconstruct a phylogeny for the region of interest
#####################################################################
## instal raxml-ng
conda install -c bioconda raxml-ng 

## run raxml-ng
raxml-ng --all --msa window12.fas --msa-format FASTA --model GTR --bs-trees 100 --threads 1 --prefix win1_and_win2




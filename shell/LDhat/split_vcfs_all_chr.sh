#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --mem=118GB

# Set the path to the LDhat directory
ldhat_dir="/.../LDhat/"

# Load required modules
module load vcftools

# Define function to perform rounding in Bash
round() {
    printf "%.${2:-0}f" "$1"
}

# Define function to perform arithmetic operations in Bash
math() {
    echo "$@" | bc -l
}

# Array of chromosome names
chroms=('Scaffold_2__2_contigs__length_236349107' 'Scaffold_1__1_contigs__length_194224858' 'Scaffold_5__2_contigs__length_177367226' 'Scaffold_3__2_contigs__length_170906183' 'Scaffold_6__2_contigs__length_157086139' 'Scaffold_7__2_contigs__length_128225229' 'Scaffold_4__1_contigs__length_105873206' 'Scaffold_10__2_contigs__length_97064778' 'Scaffold_8__1_contigs__length_80766577' 'Scaffold_9__1_contigs__length_67344957'  'Scaffold_12__2_contigs__length_57092889''Scaffold_11__1_contigs__length_50805362')

# Loop over each chromosome
for j in 0 1 2 3 4 5 6 7 8 9 10 11; do
    chrom_name=$(echo "${chroms[j]}" | awk 'BEGIN{FS="_"} {print $NF}')
    output_folder="$ldhat_dir/chr$((j+1))"
    cd "$output_folder"

    # Calculate the number of fragments
    d=$(round $(math "${chrom_name}/449000") 0)
    e=$((d + 1))

    for i in $(seq 1 $e); do
        vcftools --vcf "robeensis_highcov_chr$((j+1)).recode.vcf" \
                 --chr "${chroms[j]}" \
                 --from-bp $((449000 * (i-1) + 1)) \
                 --to-bp $((449000 * i + 1000)) \
                 --ldhat-geno \
                 --out "chr$((j+1))_highcov_$i"
    done
done



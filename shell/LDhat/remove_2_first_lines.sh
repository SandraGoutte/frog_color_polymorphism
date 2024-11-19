#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --mem 118GB


for j in {1..12}; do
    cd chr$j"_rho_pos_joint"

    for file in *_highcov_*.joint.txt; do
        gawk -i inplace 'NR>2' "$file"
    done
    cd ..
done
    
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --mem 118GB

for j in {1..12}; do
    cd chr$j"_rho_pos_joint"
    
    wc -l "chr"$j"_highcov_1.joint.txt" 

## rowbind all the tables to have a single file
    for ((i = 1; i < total_files; i++)); do
        cat "chr"$j"_highcov_"$i".joint.txt" >> "chr"$j"_highcov_1.joint.txt" 
    done
    cd ..

done    

wc -l "chr"$j"_highcov_1.joint.txt" 
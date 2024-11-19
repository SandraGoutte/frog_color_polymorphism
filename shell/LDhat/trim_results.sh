#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --mem 118GB

for j in {1..12}; do
    cd chr$j"_rho_pos_joint"

## multiply the position values by 1000 to have them in bp and not in kb
    for file in *_highcov_*.joint.txt; do
        gawk -i inplace '{$6=$6*1000;print}'  "$file"
    done  

## remove 500 bp at the end of first file and at the begining of the last file
    gawk -i inplace '$6 <= 449500' chr$j"_highcov_1.joint.txt" 

    num_files=$(ls -1 *_highcov_*.joint.txt | wc -l)
    start=($numfiles*44900)+500
    gawk -i inplace '$6 >= $start; print' chr$j"_highcov_"$num_files".joint.txt" # value is (113*44900)+500

## remove 500 bp each side of each file because they overlap by 1000 bp
    for i in {2..$num_files-1}
    do   
        min= echo $(($((449000*($i-1)))+501))
        max= echo $((( 449000*$i)+500 ))
        gawk -i inplace '{$6 >= $min; print}' "chr"$j"_highcov_"$i".joint.txt "
        gawk -i inplace '{$6 <= $max; print}' "chr"$j"_highcov_"$i".joint.txt "
    done    

## rowbind all the tables to have a single file
    for i in {2..$num_files}
    do   
        cat "chr"$j"_highcov_"$i".joint.txt" >> "chr"$j"_highcov_1.joint.txt" 
    done
done
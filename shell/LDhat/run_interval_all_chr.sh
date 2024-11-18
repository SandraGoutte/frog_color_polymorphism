#!/bin/bash
#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --time=96:00:00
#SBATCH --mem 118GB

round() {
    printf "%.${2:-0}f" "$1"
}

math() {
    echo "$*" | bc -l
}

chroms=('Scaffold_2__2_contigs__length_236349107' 'Scaffold_1__1_contigs__length_194224858' 'Scaffold_5__2_contigs__length_177367226' 'Scaffold_3__2_contigs__length_170906183' 'Scaffold_6__2_contigs__length_157086139' 'Scaffold_7__2_contigs__length_128225229' 'Scaffold_4__1_contigs__length_105873206' 'Scaffold_10__2_contigs__length_97064778' 'Scaffold_8__1_contigs__length_80766577' 'Scaffold_9__1_contigs__length_67344957' 'Scaffold_12__2_contigs__length_57092889' 'Scaffold_11__1_contigs__length_50805362' )

for j in 0 1 2 3 4 5 6 7 8 9 10 11
do 
cd /.../LDhat/chr$(($j+1))

d=$(round $(math "$(echo ${chroms[j]} | awk 'BEGIN{FS="_"} {print $NF}')/449000") 0)
e=$(($d+1))

for i in $( seq 115 $e )
do  
/.../LDhat/interval -seq chr$(($j+1))_highcov_$i.ldhat.sites -loc chr$(($j+1))_highcov_$i.ldhat.locs  -lk /.../LDhat/new_lk_20seq.txt -its 10000000 -samp 5000 -bpen 5 -prefix chr$(($j+1))_highcov_$i.
done   

done

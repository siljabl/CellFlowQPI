#!/bin/bash


# Size of DIC-box (number of neighbors on each side)
if [ $# -eq 0 ]
    then
    read -p "Enter SzW (number of neighbors on each side) : " SzW
else
    SzW=$1
fi

# Read inputs from file
while read line
do
    [[ "$line" =~ ^#.*$ ]] && continue  # Skip line starting with #
    args[$index]=$line
    index=$(($index+1))
done < config.txt

folder=${args[0]}   # path to data folder
prename=${args[1]}  # filename w.o. frame number

# Define frames
f_0=${args[2]}      # reference frame
f_MIN=$(($f_0 + 1)) # first moving frame
f_MAX=${args[3]}    # final moving frame

# Previous DIC runs to delete. Where is it specified that DIC does 6 runs?
run_MAX=5

# Change and create dictionaries
cd $folder
mkdir -p dic
mkdir -p dic/tmp_SzW_$SzW
mkdir -p dic/x_displacement_SzW_$SzW
mkdir -p dic/y_displacement_SzW_$SzW
cd tif              # must run from same folder as data

for frame in $(seq $f_MIN $f_MAX)
do
    mm3d MM2DPosSism $prename$f_0.tif $prename$frame.tif DirMEC=../dic/tmp_SzW_$SzW/$frame/ SzW=$SzW
	
    # Removing irrelevant data
    rm ../dic/tmp_SzW_$SzW/$frame/Correl_*
    rm ../dic/tmp_SzW_$SzW/$frame/Makefile*
    rm ../dic/tmp_SzW_$SzW/$frame/Masq*
    rm ../dic/tmp_SzW_$SzW/$frame/*.xml
    rm ../dic/tmp_SzW_$SzW/$frame/Z_N*
    rm ../dic/tmp_SzW_$SzW/$frame/TA*
    rm ../dic/tmp_SzW_$SzW/$frame/Px1_000*
    rm ../dic/tmp_SzW_$SzW/$frame/Px2_000*
    rm -r Pyram

    # Removing tmp from previous DIC runs
    for run in $(seq 1 $run_MAX)
    do
	rm ../dic/tmp_SzW_$SzW/$frame/Px1_Num$run*
	rm ../dic/tmp_SzW_$SzW/$frame/Px2_Num$run*
    done

    mv ../dic/tmp_SzW_$SzW/$frame/Px1_Num6_DeZoom1_LeChantier.tif ../dic/x_displacement_SzW_$SzW/$prename$f_0-$frame.tif
    mv ../dic/tmp_SzW_$SzW/$frame/Px2_Num6_DeZoom1_LeChantier.tif ../dic/y_displacement_SzW_$SzW/$prename$f_0-$frame.tif
    
done

rm -r ../dic/tmp_SzW_$SzW
mv mm3d-LogFile.txt ../dic/SzW_$SzW-LogFile.txt



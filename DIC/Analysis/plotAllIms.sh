#!/bin/bash
fields=( uniform curl field_1 field_2 )
umax=( 40 80 )

for field in "${fields[@]}"; do
for u     in "${umax[@]}";   do

python compare_ims_DIC_theoretical.py $field/umax_$u

done
done


# make this run for also frame?

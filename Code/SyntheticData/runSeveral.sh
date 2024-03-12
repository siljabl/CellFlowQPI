#!/bin/bash

fields=( uniform curl field_1 field_2 )
umax=( 40 80 )
n_frames=2
pad_width=60

for field in "${fields[@]}"; do
for u     in "${umax[@]}";   do

python generate_data.py ../../Data/InitialConditions/Well1-1_resc_reg2.tif ../../Data/SyntheticTestData/$field/ -n_frames=$n_frames -pad_width=$pad_width -u_max=$u

done
done

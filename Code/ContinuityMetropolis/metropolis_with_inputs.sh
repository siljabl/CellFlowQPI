#!/bin/bash

python tiff_to_txt.py "../../Data/SyntheticTestData/curl/umax_12/tif/Well1-1_resc_reg2_0.tif" "input/image1"
python tiff_to_txt.py "../../Data/SyntheticTestData/curl/umax_12/tif/Well1-1_resc_reg2_1.tif" "input/image2"
python tiff_to_txt.py "../../Data/SyntheticTestData/curl/umax_12/dic/x_displacement_SzW_20/Well1-1_resc_reg2_0-1.tif" "input/dx"
python tiff_to_txt.py "../../Data/SyntheticTestData/curl/umax_12/dic/y_displacement_SzW_20/Well1-1_resc_reg2_0-1.tif" "input/dy"

image_x_res=784
image_y_res=784
seed1=21;
nosteps=100000000;
alpha=1.0;
temp=70.0;
phi1=1.0;
phi2=1.0;
gamma1=1.0;
gamma2=0.0286;

./execute.x "input/image1.txt" "input/image2.txt" "input/dx.txt" "input/dy.txt" $image_x_res $image_y_res $seed1 $nosteps $alpha $temp $phi1 $phi2 $gamma1 $gamma2
python plot_output.py
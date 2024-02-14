#!/bin/bash
ss=2
prename=Well1-1_resc_reg
mkdir results_av_$ss
cd RescaledRegistered_sm



for i in {97..100}
do
    mm3d MM2DPosSism Well1-1_resc_reg2.tif $prename$i.tif DirMEC=../results_av_$ss/$i/ SzW=$ss
    
    rm ../results_av_$ss/$i/Correl_*
    rm ../results_av_$ss/$i/Makefile*
    rm ../results_av_$ss/$i/Masq*
    rm ../results_av_$ss/$i/*.xml
    rm ../results_av_$ss/$i/Z_N*
    rm ../results_av_$ss/$i/TA*
    rm ../results_av_$ss/$i/Px1_000*
    rm ../results_av_$ss/$i/Px2_000*
    rm -r Pyram
    for j in {1..5}
    do
	rm ../results_av_$ss/$i/Px1_Num$j*
	rm ../results_av_$ss/$i/Px2_Num$j*
    done
    
done


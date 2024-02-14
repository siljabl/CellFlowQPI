#!/bin/bash


name=_av_2
mkdir delta_X$name
mkdir delta_Y$name

for i in {0..100}
do
    mv results$name/$i/Px1_Num6_DeZoom1_LeChantier.tif delta_X$name/$i.tif
    mv results$name/$i/Px2_Num6_DeZoom1_LeChantier.tif delta_Y$name/$i.tif
done

rm -r stack/pyram

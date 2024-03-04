#!/bin/bash

fields=( none uniform curl field_1 field_2 )
umax=( 4 )
SzW=( 1 2 5 10 20 50 )

for field in "${fields[@]}"; do
for u     in "${umax[@]}";   do
# updating config file
sed -i -e "2s#.*#../Data/SyntheticTestData/$field/u_max_$u#" config.txt

for szw in "${SzW[@]}";
do
  #./runDIC.sh $szw
  echo $szw
done

done
done


# make this run for also frame?

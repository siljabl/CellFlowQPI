#!/bin/bash

SzW=( 1 2 5 10 20 50 )

length=${#SzW[@]}

for (( j=0; j<length; j++ ));
do
  ./runDIC.sh ${SzW[$j]}
done


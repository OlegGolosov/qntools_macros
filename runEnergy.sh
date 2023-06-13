#!/bin/bash

in=$1
p=$2

for ref in psd1 psd2 psd3 psd1_psd2_psd3 ;do
#  for comp in X Y X_Y ;do
  for comp in X ;do
    out=${in/corr.root/plots}
    out=${out}_${ref}_${comp}.root
    root -b -l -q plotFlow.C"(\"$in\",\"$out\",$p,\"$ref\",\"$comp\")"
  done 
done 

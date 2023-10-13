#!/bin/bash 

ONESAMP2COAL_MINALLELEFREQUENCY=0.05
mutationRate="0.000000012"
rangeNe=100,500
theta=0.000048,0.0048
microsatsOrSNPs=s
NeVals="00500"
numPOP="00256"

outputSampleSizes=(50 100 200)
locis=(40 80 160 320)

for outputSampleSize in "${outputSampleSizes[@]}"; do
  for loci in "${locis[@]}"; do
    for i in {1..30}; do
      ./refactor -t1 -rC -b$NeVals -d1 -u$mutationRate -v${theta} -$microsatsOrSNPs -l$loci -i$outputSampleSize -o1 -f$ONESAMP2COAL_MINALLELEFREQUENCY -p > "data/genePop${outputSampleSize}x${loci}_t${i}"
    done
  done
done


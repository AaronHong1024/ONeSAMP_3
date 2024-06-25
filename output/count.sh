#!/bin/bash

for run in {1..10}; do
  cat genePop*_$run/* | grep -o '\[\([0-9]*\.[0-9]*\)\]' | tr -d '[]' | awk '{ if ($1 < 200) print $1 }' | wc -l

done

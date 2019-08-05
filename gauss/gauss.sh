#!/bin/bash

rm stdev.dat

for i in `seq 0.001 0.001 0.01`; do

  rm -rf $i
  mkdir $i
  cd $i
  cat ../gauss.py | sed "s/XXX/${i}/g" > gauss.py
  python gauss.py
  cat stdev.tmp >> ../stdev.dat
  echo $i
  cd ..
  rm -rf $i

done

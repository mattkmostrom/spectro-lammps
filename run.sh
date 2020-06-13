#!/bin/bash

rm -rf screening
mkdir screening
cd screening
#
#for m in morse harmonic; do
##for m in morse; do
#
#	rm -rf $m
#	mkdir $m
#	                               #Starting in home directory
#  cd $m                          #Model Directory
#	for n in 1.0 1.1 1.2 1.3 1.4 1.5; do
#    mkdir $n
#    cd $n
#    for v in 0.5 1.0 1.5 2.0; do
#                                 #Temp directory
#      mkdir $v                   #Angular directory
#        cd $v
#        cp ../../../../src/*.py .
#        cp ../../../../models/${m}.in .
#        cat ../../../../models/morse.init | sed "s/XXX/$n/g" > tmp.init
#        cat tmp.init | sed "s/YYY/$v/g" > morse.init; rm tmp.init
#        python mdlammps.py ${m}.in
#        python vac.py vel.dat vac.dat
#        cd ..	                  #Back to temp
#      echo "done with ${m}_${n}_${v}"
#    done
#    echo ""
#    cd ..                       #Back to Model
#  done
#	cd ..                         #Back to home
#done
#
#a=""; for x in $(find . -name "dct.dat"); do a=$a" "$x; done; echo $a;
#b=""; for x in $(find . -name "eig.dat"); do a=$a" "$x; done; echo $a;

for m in lennard-jones; do

	rm -rf $m
	mkdir $m
	                               #Starting in home directory
  cd $m                          #Model Directory
	for n in 1.0; do
    mkdir $n
    cd $n
    for v in 0.5; do
                                 #Temp directory
      mkdir $v                   #Angular directory
        cd $v
        cp ../../../../src/*.py .
        cp ../../../../models/${m}.in .
        cat ../../../../models/morse.init | sed "s/XXX/$n/g" > tmp.init
        cat tmp.init | sed "s/YYY/$v/g" > morse.init; rm tmp.init
        python mdlammps.py ${m}.in
        python vac.py vel.dat vac.dat
        cd ..	                  #Back to temp
      echo "done with ${m}_${n}_${v}"
    done
    echo ""
    cd ..                       #Back to Model
  done
	cd ..                         #Back to home
done

a=""; for x in $(find . -name "dct.dat"); do a=$a" "$x; done; echo $a;
b=""; for x in $(find . -name "eig.dat"); do b=$b" "$x; done; echo $b;

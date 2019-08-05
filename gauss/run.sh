#!/bin/bash


#for m in morse harmonic; do
for m in harmonic; do

	rm -rf $m
	mkdir $m
	cp ../models/${m}.in $m
	cp ../models/morse.init $m
  cd $m

	python mdlammps.py ${m}.in
	python moment.py
  python vac.py moment.dat
  cd ..	
		
done

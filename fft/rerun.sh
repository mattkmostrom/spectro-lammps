#!/bin/bash

#for m in morse harmonic; do

for n in 1.2 1.3 1.4 1.5 1.6; do
 
  rm -rf $n
  mkdir $n
	cp moment.py $n
	cp ${n}.dat ${n}/fft.dat
  cd $n
	python moment.py
  cd ..
		
  echo ""
	echo "done with $n"
	echo ""	
	#cat ${n}/moment.dat | sed "s/XXX/$n/g" >> ../foot.dat
	cat ${n}/fft_moment.dat | sed "s/XXX/$n/g" >> ../foot.dat
    	
done
cd ..

echo "# distance, first, and second moments" >> head.tmp
cat foot.dat >> head.tmp
mv head.tmp fft_moment.dat

rm foot.dat

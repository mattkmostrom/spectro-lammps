#!/bin/bash

rm *.dat
#for m in morse harmonic; do
for m in morse; do


	rm -rf $m
	mkdir $m
	cp *.py $m
	cp ${m}.in $m
	cp morse.init $m
	cp moment.py $m
  cd $m	

	for n in 1.2 1.3 1.4 1.5 1.6; do

		rm -rf $n
		mkdir $n
		cat morse.init | sed "s/XXX/$n/g" > $n/morse.init
		cp *.py $n
		cp ${m}.in $n
		cd $n
		python mdlammps.py ${m}.in
		python vac.py vel.dat vac.dat
		python moment.py
		#mv moment.dat ${m}_${n}_moment.dat
    cd ..	
		
    echo ""
		echo "done with $n"
		echo ""	
		cat ${n}/moment.dat | sed "s/XXX/$n/g" >> ../foot.dat
    	
  done
	cd ..
done

echo "# distance, first, and second moments" >> head.tmp
cat foot.dat >> head.tmp
mv head.tmp moment.dat

a=""; for x in $(ls */*/vac.dat); do a=$a" "$x; done; 
a=""; for x in $(ls */*/*eig.dat); do a=$a" "$x; done;  
xmgrace $a

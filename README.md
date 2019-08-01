# spectro-lammps

This is a code written primarily by Dr. Preston Moore, Andres Saez, and Matthew
Mostrom in Python for the purpose of using classical dynamics to characterize
spectroscopic phenomena of molecular systems. Input currently is made to
resemble LAMMPS input scripts for the purpose of being able to quickly and
easily check the validity of our results with the folks at Sandia for any system
we care to look at.

Running MD:
Like LAMMPS, our code takes the name of the program to be run, followed by input
file. A sample input line would look like this:

python mdlammps.py morse.in

#!/usr/bin/python3

import math

# loop over grid and plot potentials

gridsize = 50
r0 = 1.5
k = 2.2
xmin = -r0*.25
xmax = r0*1.5
dx = (xmax-xmin)/(gridsize-1)
for i in range(gridsize):
    x = i*dx+xmin
    for j in range(gridsize):
        y = j*dx+xmin
        r = math.sqrt(x**2+y**2)
        v = k*(r-r0)**2
        print(x,y,v)
    print()
#print()

exit(0)
D = 1.5
alpha = 0.75
r0 = 1.5
xmin = r0*.25
xmax = r0*2.5
dx = (xmax-xmin)/(gridsize-1)

for i in range(gridsize):
    x = i*dx+xmin
    for j in range(gridsize):
        y = j*dx+xmin
        r = math.sqrt(x**2+y**2)
        dr = r-r0
        expar = math.exp(-alpha*dr)
        exparm1 = 1-expar
        pot = D*exparm1*exparm1
        print(x,y,pot)
    print()

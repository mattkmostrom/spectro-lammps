#!/usr/bin/python3

import numpy
import math
from scipy.fftpack import dct

x0 = 0.
A = 1.
fit = XXX # Alpha; Arbitrary fitting parameter
dt = 0.2

compute = numpy.zeros(3000)
gauss = numpy.zeros(len(compute))
nwindow = len(compute)

for i in range(len(compute)):
    gauss[i] = float(A*math.exp(-fit*((0.5*dt*i-x0)**2))) 


#gauss = []
#for i in range(len(compute)):
#    compute = float(A*math.exp(-fit*((0.5*dt*i-x0)**2))) 
#    print(compute)
#    gauss.append(compute.real)

#print("Writing gauss.dat")
file = open("gauss.dat","w")
for i in range(len(gauss)):
    line = str(dt*i) + " " + str(gauss[i]) + "\n"
    file.write(line)

gauss_dct = dct(gauss) #Cosine transform
gauss_dct = gauss_dct*(1/gauss_dct[0]) #Normalization

#print("Writing gauss_dct.dat")
file = open("gauss_dct.dat","w")
for i in range(len(gauss)):
    line = str(i*math.pi/(dt*nwindow)) + " " + str(gauss_dct[i]) + "\n"
    file.write(line)

rangemax = 0.65
rangemin = 0.6
file = open("gauss_dct.dat","r")
lines = file.readlines()

stdev = []
for i in lines:
    if rangemax > float(i.split()[1]) > rangemin:
        stdev.append(i)

line = str("XXX")+" "+str(stdev)+"\n"

print("St. Dev Domain: ",stdev)

file = open("stdev.tmp","w")
file.write(line)

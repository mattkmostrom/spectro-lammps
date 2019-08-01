#takes a file of two columns, x,y sorted by x without comments at the top

#need a way to specify which x-values to search
#need a way to specify which moment.dat is from which distance
from __future__ import print_function
import numpy as np
from scipy.integrate import simps
from numpy import trapz
import glob

file = open("morse.eig.dat","r") # Figure out how to use wildcard

lines = file.readlines()
total = np.zeros((len(lines)))
data = np.zeros((len(lines)))
print("")

for i in range(1,len(lines)): # Get all of the x-values
    total[i] = (lines[i]).split()[0]

domain = []
limit = 1.0
for i in total: # Filter out everything left of "limit"
    if i and i > limit:
        domain.append(i)

for i in range(len(domain)): # Put y's associated with domain into data 
    data[i] = float((lines[-(i+1)]).split()[1])

threshold = 0.05
samples = []
for i in data: # Filter out signal below "threshold"
    if i and i > threshold:
        samples.append(i)

print("Data Part 2: ",samples)
print("")

dt = total[2]-total[1] # Calculate integration step; changes with bin count
print("dx: ",dt)
xmin = domain[0] #First one
print("xmin: ",xmin)
xmax = domain[-1] # Last one
print("xmax: ",xmax)

# composite trapezoidal rule.
area = trapz(samples, dx=dt)
mom1 = (1/(xmax-xmin)*area)
print("First Moment: ", mom1)

mom2 = np.std(data)
print("Second Moment: ", mom2)

file = open("moment.dat","w")
line = str("XXX")+" "+str(mom1)+" "+str(mom2)+ "\n"
file.write(line)

#sample data
#1.38149 0.0869713
#1.40832 0.186367
#1.43515 0.173943
#1.46198 0.198791

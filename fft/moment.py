#takes a file of two columns, x,y sorted by x without comments at the top
from __future__ import print_function
import math
import numpy as np
from scipy.integrate import simps
from numpy import trapz

file = open("fft.dat","r") # Figure out how to use wildcard

lines = file.readlines()
total = np.zeros((len(lines)))
print("")

for i in range(len(lines)): # Get all of the x-values
    total[i] = (lines[i]).split()[0]

domainmin = 0.      # Left-most value
domainmax = 0.075   # Right-most value
domain = []
for i in total: # Filter out everything left of "limit"
    if i and domainmax > i > domainmin:
        domain.append(i)

data = np.zeros((len(domain)))

for i in range(len(domain)): # Put y's associated with domain into data 
    #data[i] = float((domain[-i]*lines[-(i+1)]).split()[1])
    data[i] = float((lines[(i)]).split()[1]) # Can only use if there is no header

print("")
print("total: ", total)
print("")
print("Domain: ",domain)
print("")
print("Data: ",data)

rangemin = 0.0      #signal threshold
sig0 = []   # f_x
sig1 = []   # x*f_x
sig2 = []   # x^2*f_x
for i in range(len(data)): # Filter out signal below "threshold"
    if data[i] and data[i] > rangemin:
        sig0.append(data[i])
        sig1.append(domain[i]*data[i])
        sig2.append(domain[i]*domain[i]*data[i])

print("")
print("Samples: ",sig1)
print("")

dt = total[2]-total[1] # Calculate integration step; changes with bin count
print("dx: ",dt)
xmin = domain[0] #First one
print("xmin: ",xmin)
xmax = domain[-1] # Last one
print("xmax: ",xmax)

# composite trapezoidal rule.
area = trapz(sig1, dx=dt)
integral = trapz(sig0, dx=dt)
extra_int = trapz(sig2, dx=dt)

print("Distribution Integral: ", integral)

mom1 = ((1/(integral))*area)
print("First Moment: ", mom1)

mom2 = 2*math.sqrt(((1/(integral))*extra_int) - mom1**2) # Variance
print("Second Moment: ", mom2)

file = open("fft_moment.dat","w")
line = str("XXX")+" "+str(mom1)+" "+str(mom2)+ "\n"
file.write(line)

#sample data
#1.38149 0.0869713
#1.40832 0.186367
#1.43515 0.173943
#1.46198 0.198791

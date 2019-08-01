#!/usr/bin/python2
# run MD using lammps style input

import sys
import numpy as np
import math
import time
import re

#import mdglobal  # file with global variables
import mdinput   # file with input routines
import mdoutput  # file with output routines
import mdbond    # file with bonding routines
import mdstep    # file with the integration

# assigned global variables
global natoms       # number of atoms
global atypes       # number of atom types
global nbonds       # number of bonds
global tbonds       # number of bond types
global box          # box coordinates
global pot          # potential components
global mass         # masses of each type
global masses       # broadcast of mass type to atoms
global pos          # positions of each atom
global vel          # velocities of each atom
global acc          # acceleration of each atom
global aatype       # array of atom types
global bonds        # bonds (array with type, ibond, jbond)
global hessian      # hessian matrix
global abtype       # array of bond types
global logfile      # file to output thermodata
global fix_type

box = np.zeros(3)
pot = np.zeros(6)
fix_type = None

kb = 1.38064852e-23
zeta = np.zeros(2)
vtherm = np.zeros(2)
G = np.zeros(2)

#------------------------------------------------

def readin(): # read lammps like infile

    global nsteps, dt, initfile, ithermo, idump, dumpfile, bond_style, bondcoeff
    global logfile, inmfile, inmo
    global bondcoeff, reps, fix_type
    global natoms, atypes, nbonds, tbonds, box

    # print lines
    data, bond_style, bondcoeff, reps, fix_type, var_lst = mdinput.readsysvals(sys.argv[1]) # read lammps in file
    dt, initfile, bond_styles, idump, dumpfile, ithermo, logfile, inmfile, inmo, nsteps = data
    print("dt, initfile, bond_styles, idump, dumpfile, ithermo, logfile, inmfile, inmo, nsteps",data)
    
    natoms, atypes, nbonds, tbonds, box[0], box[1], box[2] = mdinput.readinvals(initfile) 
    print("Natoms",natoms," Atypes",atypes," Bonds",nbonds," Btypes",tbonds)
    print("Box",box)

    # allocate arrays from data
    global mass, aatype, pos, vel, acc, masses, bonds, hessian, zeta, Q

    acc = np.zeros((natoms,3))

    mass, aatype, pos, vel, masses, bonds = mdinput.make_arrays(initfile,reps)
    
    if re.search('nvt',fix_type,flags=re.IGNORECASE):
        global T, Tdamp, Q
        var = [float(num) for num in var_lst[1:4]]
        T, T, Tdamp = var
        Q = np.array([3*natoms*T*Tdamp*Tdamp]*2)

#-----------------------------------------------------------
def force(pos,vel,acc,masses,nbond,bonds,bondcoeff,pot): # get forces from potentials

    acc.fill(0) # zero out forces/acceration
    # lj
    pot[0] = 0
    # bonds
    pot[1] = mdbond.bond_force(bond_style,nbonds,bonds,bondcoeff,pos,acc,masses)
    # bend
    pot[2] = 0
    # torsion
    pot[3] = 0

    # Change forces into accelerations
    acc /= masses

    return
#-----------------------------------------------------------

# read command line for input file
if (len(sys.argv) < 2):  # error check that we have an input file
    print("No input file? or wrong number of arguments")
    exit(1)

if not re.search(re.compile(r'.+\.in'),sys.argv[1]):
    print('Incorrect input file type.')
    exit(1)

print (sys.argv)

readin() # read infile

ke = (0.5*np.dot(masses.transpose()[0],np.array([np.dot(vec,vec) for vec in vel])))
ke_init = ke

if fix_type == 'nvt':
    def step(): # nose-hoover chain
        global pos, vel, acc, dt, ke, kb, w
        
        ke,vel = mdstep.nhchain(Q,G,dt,natoms,vtherm,zeta,ke,vel,T)
        vel += acc*dt/2.0
        pos += vel*dt
        force(pos,vel,acc,masses,nbonds,bonds,bondcoeff,pot) 
        vel += acc*dt/2.0
        ke,vel = mdstep.nhchain(Q,G,dt,natoms,vtherm,zeta,ke,vel,T)
elif fix_type == 'nve':
    def step(): # velocity verlet
        global pos, vel, acc, dt
        
        vel += acc*dt/2.0
        pos += vel*dt
        force(pos,vel,acc,masses,nbonds,bonds,bondcoeff,pot)
        vel += acc*dt/2.0    
else:
    print('Unrecognized fix type')
    exit(1)

# inital force and adjustments
mdstep.zero_momentum(masses,vel)  # zero the momentum
force(pos,vel,acc,masses,nbonds,bonds,bondcoeff,pot) 
teng = mdoutput.write_thermo(logfile,0,natoms,masses,pos,vel,pot)

itime = 1
tnow = time.time()
ttime = tnow
tol = 1e-8
dump_vel = 5

print("Running dynamics")

eig_array = [] # empty array for the eigenvalues

for istep in range(1,nsteps+1):

    step() # take a step

    if(istep%inmo==0): # get instantaneous normal modes
        hessian = mdbond.inm(bond_style,nbonds,bonds,bondcoeff,pos,masses)

        # print(hessian)
        w,v = np.linalg.eig(hessian)
        # remove 3 lowest eigegvalues (translations of entire system hopefully)
        for i in range(3): # Be careful and make sure the eigenvector is translation
            idx = np.argmin(np.abs(w.real))
            if(abs(w[idx]) > tol):
                print("Removing eigvalue > tol",w[idx])
            w = np.delete(w,idx)
        eig_array.append(w.real) # only get real part of array - imag do to round off error is small so we throw away.

    if(istep%ithermo==0): # write out thermodynamic data
        teng = mdoutput.write_thermo(logfile,istep,natoms,masses,pos,vel,pot)

    if(istep%idump==0): # dump to xyz file so we can see this in lammps
        mdoutput.write_dump(dumpfile,istep*dt,natoms,pos,aatype)

    if(istep%dump_vel==0): # dump out vel file to correlate
        mdoutput.write_dump_vel("vel.dat",istep*dt,natoms,vel)

    if(itime < time.time()-tnow): # report where we are
        print('step = {}/{} = {:.4f}%, teng = {:g}, time = {:g}'.format(istep,nsteps,istep/nsteps*100,teng,time.time()-ttime))
        tnow = time.time()

print('Done dynamics! total time = {:g} seconds'.format(time.time()-ttime))
mdoutput.write_init("test.init",istep-1,natoms,atypes,nbonds,tbonds,box,mass,pos,vel,bonds,aatype)

if(fix_type=='nve'):
    print("NVE Energy Drift = ",ke-ke_init)
elif(fix_type=='nvt'):
    therm_energy = (0.5*np.dot(Q,np.array([np.dot(vec,vec) for vec in vtherm])))
    print("NVT Energy Diff = ",ke-ke_init+therm_energy)
else:
    print('Fix_type unknown')
    
eig_array = np.array(eig_array)
eig_array = [np.sign(x)*math.sqrt(abs(x)) for x in eig_array.ravel()]

#Create histogram!
nconf = len(eig_array)
if(nconf==0):
    print("No configurations calculated eigenvalues! thus NOT calculating historgram")
else:
    print("Creating Histogram with",len(eig_array),"configurations")
    #q1, q3 = np.percentile(np.array(eig_array), [25, 75])
    #iqr = q3 - q1

    #fd_width = 2*iqr/(nconf**(1/3))
    #fd = (np.amax(eig_array) - np.amin(eig_array))/fd_width
    #fd = int(fd) + 1

    #sturges = np.log2(nconf) + 1
    #sturges = int(sturges) + 1

    #bin_ct = max(fd,sturges)
    bin_ct = 100
    histo,histedge = np.histogram(np.array(eig_array),bins=bin_ct,density=True)
    histdat = np.zeros((histo.size,2))
    for i in range(histo.size):
        histdat[i][0] = (histedge[i]+histedge[i+1])/2
        histdat[i][1] = histo[i]
        #print(histo,histedge,histdat)
    head = "Histogram of eigenvalues " + sys.argv[0] + " " + str(len(eig_array))
    np.savetxt(inmfile,(histdat),header=head,fmt="%g")

print("Done!")
exit(0)

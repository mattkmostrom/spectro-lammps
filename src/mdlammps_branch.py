#!/usr/bin/python3
# run MD using lammps style input

import sys
import numpy
import time

#import mdglobal  # file with global variables
import mdinput   # file with input routines
import mdoutput  # file with output routines
import mdlj      # file with non-bonded routines
import mdbond    # file with bonding routines

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

box = numpy.zeros(3)
pot = numpy.zeros(6)

#-------------------------------------------------
def readinit(datafile): # read lammps init data file

    global natoms, atypes, nbonds, tbonds, box

    print ("Reading",datafile)
    fi = open(datafile,"r")
    lines = fi.readlines() # read in lines all at once
    fi.close()

    data =[0]*7
    mdinput.readinvals(lines,data)

    # unpack or destructuring
    natoms, atypes, nbonds, tbonds, box[0], box[1], box[2] = data
    print("Natoms",natoms," Atypes",atypes," Bonds",nbonds," Btypes",tbonds)
    print("Box",box)

    # allocate arrays from data
    global mass, aatype, pos, vel, acc, masses, bonds, hessian

    mass = numpy.zeros(atypes)
    aatype = numpy.zeros(natoms,dtype=int)
    pos = numpy.zeros((natoms,3))
    vel = numpy.zeros((natoms,3))
    acc = numpy.zeros((natoms,3))
    masses = numpy.zeros((natoms,3))
    bonds = numpy.zeros((nbonds,3),dtype=int)
    hessian = numpy.zeros((3*natoms,3*natoms))

    mdinput.getmasses(lines,atypes,mass)
    mdinput.getatoms(lines,natoms,aatype,pos,mass,masses)
    mdinput.getvel(lines,natoms,vel)
    mdinput.getbonds(lines,nbonds,bonds)

#-------------------------------------------
def readin(): # read lammps like infile

    global nsteps, dt, initfile, ithermo, idump, dumpfile, bond_style, bondcoeff
    global logfile, inmfile, inmo

    logfile = None
    inmfile = None
    print ("Reading",sys.argv[1])
    fi = open(sys.argv[1],"r")
    lines = fi.readlines() # read in lines all at once
    fi.close()

    # print lines
    data =[None]*10
    mdinput.readsysvals(lines,data) # read lammps in file
    nsteps, dt, initfile, ithermo, idump, dumpfile, bond_styles, logfile, inmfile, inmo = data
    print("nsteps, dt, initfile, ithermo, idump, dumpfile, bond_styles, logfile, inmfile, inm",data)

    readinit(initfile)  # read initfile to get atoms bonds and types

    global bondcoeff
    bondcoeff = numpy.zeros((tbonds,3))  # allocate bondcoeff
    bond_style = mdinput.getbondcoeff(lines,bond_styles,tbonds,bondcoeff) # readin bondcoeff

#-----------------------------------------------------------
def force(): # get forces from potentials
    global pot, nbonds, bonds, bondcoeff
    global masses, pos, vel, acc

    acc.fill(0) # zero out forces/acceration
    # lj
    pot[0] = 0
    # bonds
    pot[1] = mdbond.bond(bond_style,nbonds,bonds,bondcoeff,pos,acc,masses)
    # bend
    pot[2] = 0
    # torsion
    pot[3] = 0

    # change forces into accelerations
    acc /= masses

#-----------------------------------------------------------
def step(): # velocity verlet (using 1/2 steps)
    global pos, vel, acc, dt
    # print istep,pos,vel,acc
    vel += acc*dt/2.0
    pos += vel*dt
    force()
    vel += acc*dt/2.0

#-----------------------------------------------------------

# read command line for input file
if (len(sys.argv) != 2):  # error check that we have an input file
    print("No input file? or wrong number of arguments")
    exit(1)
print (sys.argv)

readin() # read infile

# inital force and adjustments
mdlj.zero_momentum(masses,vel)  # zero the momentum
force()
teng = mdoutput.write_thermo(logfile,0,natoms,masses,pos,vel,pot)

itime = 10 # report total energy every 1 seconds
tnow = time.time()
ttime = tnow
tol = 1e-8
#inmo = 100 # write hessian ever 10 steps

print("Running dynamics")

eig_array = []

for istep in range(1,nsteps+1):

    step() # take a step

    if(istep%inmo==0): # get instantaneous normal modes
        hessian = numpy.zeros((pos.size,pos.size))
        mdbond.inm(bond_style,nbonds,bonds,bondcoeff,pos,masses,hessian)

        # print(hessian)
        w,v = numpy.linalg.eig(hessian)
        # remove lowest eigegvalues (translations of entire system)
        idx = numpy.argmin(numpy.abs(w.real))
        if(abs(w[idx]) > tol):
            print("Warning! Removing eigenvalue > tol",w[idx])
        w = numpy.delete(w,idx)
        idx = numpy.argmin(numpy.abs(w.real))
        if(abs(w[idx]) > tol):
            print("Warning! Removing eigenvalue > tol",w[idx])
        w = numpy.delete(w,idx)
        idx = numpy.argmin(numpy.abs(w.real))
        if(abs(w[idx]) > tol):
            print("Warning! Removing eigenvalue > tol",w[idx])
        w = numpy.delete(w,idx)
        eig_array.append(w.real) # only get real part of array - imag do to round off error is small so we throw away.
        # mdoutput.write_inm(istep,hessian)

    if(istep%ithermo==0): # write out thermodynamic data
        teng = mdoutput.write_thermo(logfile,istep,natoms,masses,pos,vel,pot)

    if(istep%idump==0): # dump to xyz file so we can see this in lammps
        mdoutput.write_dump(dumpfile,istep,natoms,pos,aatype)

    if(itime < time.time()-tnow): # report where we are
        print('step = {}/{} = {:.4f}%, teng = {:g}, time = {:g}'.format(istep,nsteps,istep/nsteps*100,teng,time.time()-ttime))
        tnow = time.time()

print('Done dynamics! total time = {:g} seconds'.format(time.time()-ttime))
mdoutput.write_init("test.init",istep-1,natoms,atypes,nbonds,tbonds,box,mass,pos,vel,bonds,aatype)

#Create histogram!
nconf = len(eig_array)
if(nconf==0):
    print("No configurations calculated eigenvalues! thus NOT calculating historgram")
else:
    print("Creating Histogram with",len(eig_array),"configurations")
    histo,histedge = numpy.histogram(numpy.array(eig_array),bins='auto',density=True)
    histdat = numpy.zeros((histo.size,2))
    for i in range(histo.size):
        histdat[i][0] = (histedge[i]+histedge[i+1])/2
        histdat[i][1] = histo[i]
        #print(histo,histedge,histdat)
    head = "Histogram of eigenvalues " + sys.argv[0] + " " + str(len(eig_array))
    numpy.savetxt(inmfile,(histdat),header=head,fmt="%g")

print("Done!")
exit(0)

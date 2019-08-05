
import numpy
import math

import mdbond
#--------------------# get forces from potentials
def force(natoms,pot,nbonds,bond_style,bondcoeff,bonds,masses,pos,vel,acc): 

    acc.fill(0) # zero out forces/acceration
    # lj
    pot[0] = 0
    # bonds
    pot[1] = mdbond.bond_force(bond_style,nbonds,bonds,bondcoeff,pos,acc)
    # bend
    pot[2] = 0
    # torsion
    pot[3] = 0

    # change forces into accelerations
    acc /= masses

# ----------------------------------------------------------
def zero_momentum(masses,vel):  # zero the liniar momentum
    mom = masses*vel # get momentum
    tmom = numpy.sum(mom,axis=0)/numpy.sum(masses,axis=0) # total mom/mass
    vel -= tmom # zero out

#-------------------------------------------------------
# constansts for NHC integration

#kb = 1 # for internal/lj units
kb = 1.38064852e-23
w = [1.0/(2.0 - 2.0**(1/3)),0,0]
w[2] = w[0]
w[1] = 1.0 - 2.0*w[0]

def nhchain(Q, G, dt, natoms, vtherm, zeta, ke, vel, T):
    
    M = len(zeta)  # chain length
    scale = 1.0
    for i in range(3):
        ts = w[i] * dt
        for j in range(1, M - 2):
            G[M - j] = (Q[M - j - 1] * vtherm[M - j - 1] * vtherm[M - j - 1] - kb * T) / Q[M - j - 1]
            vtherm[M - j] += G[M - j] * ts / 4.0
            vtherm[M - j - 1] *= math.exp(-vtherm[M - j] * ts / 8.0)
        vtherm[0] *= math.exp(-vtherm[1] * ts / 8.0)
        G[0] = (ke - (3.0 * natoms) * kb * T) / Q[M - 2]
        vtherm[0] += G[0] * ts / 4.0
        vtherm[0] *= math.exp(-vtherm[1] * ts / 8.0)
        scale *= math.exp(-vtherm[0] * ts / 2.0)
        ke *= math.exp(-vtherm[0] * ts)
        zeta += vtherm * ts / 2.0
        vtherm[0] *= math.exp(-vtherm[1] * ts / 8.0)
        G[0] = (ke - 3.0 * natoms * kb * T) / Q[M - 2]
        vtherm[0] += G[0] * ts / 4.0
        vtherm[0] *= math.exp(-vtherm[1] * ts / 8.0)
        for j in range(1, M - 2):
            vtherm[j] *= math.exp(-vtherm[j + 1] * ts / 8.0)
            G[j] = (Q[j - 1] * vtherm[j - 1] * vtherm[j - 1] - kb * T) / Q[j]
            vtherm[j] += G[j] * ts / 4.0
            vtherm[j] *= math.exp(-vtherm[j + 1] * ts / 8.0)
        G[M - 1] = (Q[M - 2] * vtherm[M - 2] * vtherm[M - 2] - kb * T) / Q[M - 1]
        vtherm[M - 1] += G[M - 1] * ts / 4.0

    vel *= scale
    return ke, vel

#-----------------------------------------------------------
def step(natoms,ensamble,dt,pot,nbonds,bond_style,bondcoeff,bonds,masses,pos,vel,acc): 
    # print istep,pos,vel,acc

# velocity verlet (using 1/2 steps)
    if(ensamble==1): # nvt
        nhc()
    vel += acc*dt/2.0
    pos += vel*dt
    force(natoms,pot,nbonds,bond_style,bondcoeff,bonds,masses,pos,vel,acc)
    vel += acc*dt/2.0

    if(ensamble==1): # nvt
        nhc()
    

# mdbond.py  # bonding routines for mdpython
#TODO
#Add angle potential
#Add LJ/bends/torsions

import numpy as np
import math

kb = 1.38064852e-23
w = [1.0/(2.0 - 2.0**(1/3)),0,0]
w[2] = w[0]
w[1] = 1.0 - 2.0*w[0]

def bond_force(bond_style,nbonds,bonds,bondcoeff,pos,acc,masses):

    pbond = 0
    for i in range(nbonds): 
        itype = bonds[i][0]
        posi = pos[bonds[i][1]]
        posj = pos[bonds[i][2]]

        rv = posj - posi
        r =  math.sqrt(np.dot(rv,rv))
        
        if bond_style == 0:
            k0 = bondcoeff[itype][0]
            r0 = bondcoeff[itype][1]
            dr = r - r0
            pot = k0*dr*dr
            dudr = 2*k0*(r-r0)
            F = (dudr/r)*rv

        elif bond_style == 1: #Morse
            D = bondcoeff[bonds[i][0]][0]
            alpha = bondcoeff[bonds[i][0]][1]
            r0 = bondcoeff[bonds[i][0]][2]
            dr = r-r0
            expar = math.exp(-alpha*dr)
            pot = D*(1.0 - expar)*(1.0 - expar)
            dudr = 2.0*D*alpha*expar*(1.0-expar)
            F = (dudr/r)*rv

        else:
            print ("Error in bond_style? in routine bond_force\n")
            exit(1)

        pbond += pot             # sum bond potential

        acc[bonds[i][1]] += F  # add forces to particles
        acc[bonds[i][2]] -= F
    
    return pbond

def inm(bond_style,nbonds,bonds,bondcoeff,pos,masses):

    hessian = np.zeros((pos.size,pos.size))

    for i in range(nbonds):  
        itype = bonds[i][0] 
        idx = bonds[i][1] 
        jdx = bonds[i][2]
        posi = pos[idx]
        posj = pos[jdx]

        rv = posj-posi
        r = math.sqrt(np.dot(rv,rv))

        ii = idx*3
        jj = jdx*3

        #d^2 /dx^2 U(r(x,y)) = r" U' + r'^2 U"
        #d^2/ dx dy U(r(x,y)) = d^2/dxdy r dU/dr + dr/dx dr/dy d^2 U/dr^2
        if(bond_style==0):  # Harmonic
            k0 = bondcoeff[itype][0] 
            r0 = bondcoeff[itype][1]
            dudr = 2*k0*(r-r0)
            du2dr2 = 2*k0

        if(bond_style==1): #Morse
            D = bondcoeff[bonds[i][0]][0]
            alpha = bondcoeff[bonds[i][0]][1]
            r0 = bondcoeff[bonds[i][0]][2]
            dr = r-r0
            expar = math.exp(-alpha*dr)
            dudr = 2.0*D*alpha*expar*(1.0-expar)
            du2dr2 = (2.0*D*alpha*alpha)*(2*expar*expar - expar)

        for k in range(3): #populate upper half of hessian 
            diagelm = dudr/r
            hessian[ii+k][ii+k] += diagelm
            hessian[ii+k][jj+k] -= diagelm
            hessian[jj+k][jj+k] += diagelm
            for l in range(3): 
                elmij = -(rv[k]*rv[l])*dudr/(r**3) + du2dr2*rv[k]*rv[l]/(r**2)
                hessian[ii+k][ii+l] += elmij
                hessian[ii+k][jj+l] -= elmij
                hessian[jj+k][jj+l] += elmij

    # mass weight
    ma = masses.reshape(pos.size)
    for i in range(pos.size):
        hessian[i][i] /= ma[i]
        for j in range(i+1,pos.size):
            hessian[i][j] /= math.sqrt(ma[i]*ma[j])
            hessian[j][i] = hessian[i][j]

    return hessian

#-------------------------------------------------------------------------

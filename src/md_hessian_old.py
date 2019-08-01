#        x0 = ipos[0]
#        y0 = ipos[1]
#        z0 = ipos[2]
#        x1 = jpos[0]
#        y1 = jpos[1]
#        z1 = jpos[2]
         #r = math.sqrt((x1-x0)**2+(y1-y0)**2+(z1-z0)**2)

#        drdx = rv[0]/r
#        drdy = rv[1]/r
#        drdz = rv[2]/r
#        dr2dxz = -(rv[0]*rv[2])/r**3
#        dr2dyz = -(rv[1]*rv[2])/r**3
#        dr2dxy = -(rv[0]*rv[1])/r**3
#        dr2dx2 = (1/r)-(rv[0]*rv[0])/r**3
#        dr2dy2 = (1/r)-(rv[1]*rv[1])/r**3
#        dr2dz2 = (1/r)-(rv[2]*rv[2])/r**3

#        hessian[0][0] = dr2dx2*dudr + drdx*drdx*du2dr2
#        hessian[0][1] = dr2dxy*dudr + drdx*drdy*du2dr2
#        hessian[0][2] = dr2dxz*dudr + drdx*drdz*du2dr2
#        hessian[0][3] = -dr2dx2*dudr - drdx*drdx*du2dr2
#        hessian[0][4] = -dr2dxy*dudr - drdy*drdx*du2dr2
#        hessian[0][5] = -dr2dxz*dudr - drdz*drdx*du2dr2

#        hessian[1][1] = dr2dy2*dudr + drdy*drdy*du2dr2
#        hessian[1][2] = dr2dyz*dudr + drdy*drdz*du2dr2
#        hessian[1][3] = -dr2dxy*dudr - drdx*drdy*du2dr2
#        hessian[1][4] = -dr2dy2*dudr - drdy*drdy*du2dr2
#        hessian[1][5] = -dr2dyz*dudr - drdz*drdy*du2dr2

#        hessian[2][2] = dr2dz2*dudr + drdz*drdz*du2dr2
#        hessian[2][3] = -dr2dxz*dudr - drdx*drdz*du2dr2
#        hessian[2][4] = -dr2dyz*dudr - drdy*drdz*du2dr2
#        hessian[2][5] = -dr2dz2*dudr - drdz*drdz*du2dr2

#        hessian[3][3] = dr2dx2*dudr + drdx*drdx*du2dr2
#        hessian[3][4] = dr2dxy*dudr + drdy*drdx*du2dr2
#        hessian[3][5] = dr2dxz*dudr + drdz*drdx*du2dr2

#        hessian[4][4] = dr2dy2*dudr + drdy*drdy*du2dr2
#        hessian[4][5] = dr2dyz*dudr + drdy*drdz*du2dr2

#        hessian[5][5] = dr2dz2*dudr + drdz*drdz*du2dr2

#    print(hessian)
def bond_harm(nbonds,bonds,bondcoeff,pos,acc):
    #Harmonic	k*(r-r0)^2
    #d/dr	2*k*(r-r0)
    #d^2/dr^2   2*k

    pbond = 0
    for i in range(nbonds):  # loop over bonds,
        ipos = pos[bonds[i][1]]
        jpos = pos[bonds[i][2]]
        bondk = bondcoeff[bonds[i][0]][0] # use type to bond params
        r0 = bondcoeff[bonds[i][0]][1]

        dpos = jpos-ipos
        r =  math.sqrt(numpy.dot(dpos,dpos))
        dr = r-r0

        pot = bondk*dr**2
        pbond += pot             # total bond potential

        dudr = 2.*bondk*dr
        dpos = (dudr/r)*dpos
        # du2dr2 = 2.*bondk

        acc[bonds[i][1]] += dpos  # add forces back
        acc[bonds[i][2]] -= dpos

    return(pbond)

#----------------Morse potential---------------------------
def bond_morse(nbonds,bonds,bondcoeff,pos,acc):
    # Morse EQ	D*(1-exp(-a(r-r0)))^2
    # d/dr	2 D a exp(-a(r-r0))*(1-exp(-a(r-r0)))
    # d^2/dr^2	4 a^2 D exp(-2 a (r-r0)) - 2 a^2 D exp(-a (r-r0))

    pbond = 0
    for i in range(nbonds):  # loop over bonds,
        ipos = pos[bonds[i][1]]
        jpos = pos[bonds[i][2]]
        # use type to bond params
        D = bondcoeff[bonds[i][0]][0]
        alpha = bondcoeff[bonds[i][0]][1]
        r0 = bondcoeff[bonds[i][0]][2]

        dpos = jpos-ipos
        r =  math.sqrt(numpy.dot(dpos,dpos))
        dr = r-r0

        expar = math.exp(-alpha*dr)
        pot = D*(1-expar)**2
        pbond += pot             # total bond

        dudr = 2.0*D * alpha * expar *(1.0-expar)
        dpos = (dudr/r)*dpos

        acc[bonds[i][1]] += dpos  # add forces back
        acc[bonds[i][2]] -= dpos
    return(pbond)


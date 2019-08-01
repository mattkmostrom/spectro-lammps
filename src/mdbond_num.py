# Routines to check the forces and hessian numerically

#--------------Numerical forces  dU/dx -------------------------------
def bond_num(bond_style,nbonds,bonds,bondcoeff,pos,acc):

    pbond = 0

#    print("Analytical forces",pos.size)
#    print(acc)
    print("Calculating Numerical Forces")
    acc_num = numpy.zeros(pos.size)
    acc_d = numpy.zeros(acc.shape)
    h = .000001

    post = pos.reshape(pos.size)
    # print(post.reshape((-1,3)))
    # uses df(x)/dx = (f(x+h)-f(x-h))/(2h)
    pbond = bond_force(bond_style,nbonds,bonds,bondcoeff,post.reshape(-1,3),acc_d)
    for i in range(post.size):
        tpos = post[i]
        post[i] = tpos + h
        pbond1 = bond_force(bond_style,nbonds,bonds,bondcoeff,post.reshape(-1,3),acc_d)
        post[i] = tpos - h
        pbondm1 = bond_force(bond_style,nbonds,bonds,bondcoeff,post.reshape(-1,3),acc_d)
        post[i] = tpos

        acc_num[i] = -(pbond1-pbondm1)/(2.0*h)
        #print(i,pbond,pbond1,pbondm1)

    numpy.copyto(acc,acc_num.reshape(-1,3))
    return pbond

#-----------------------Numerical second derivative --------------------
def bond_hess_num(bond_style,nbonds,bonds,bondcoeff,pos,acc,masses,hess_num):

    pbond = 0
    h = .0001
    ma = masses.reshape(pos.size)
    acc_d = numpy.zeros(acc.shape)

    print("Calculating Numerical 2nd derivatives")
    post = pos.reshape(pos.size)
    # print(post.reshape((-1,3)))
    for i in range(post.size):
        tpos = post[i]
        # uses d^2f(x)/ dx dx = (f(x+h,y)+f(x-h,y)-2f(x,y))/(hh)
        pbond = bond_force(bond_style,nbonds,bonds,bondcoeff,post.reshape(-1,3),acc_d)
        post[i] = tpos - h
        pbondm1 = bond_force(bond_style,nbonds,bonds,bondcoeff,post.reshape(-1,3),acc_d)
        post[i] = tpos + h
        pbond1 = bond_force(bond_style,nbonds,bonds,bondcoeff,post.reshape(-1,3),acc_d)
        post[i] = tpos

        hess_num[i][i] = (pbond1+pbondm1-2.0*pbond)/(h*h*ma[i]) #diagonal with mass weight

        #print(i,pbond1,pbondm1,pbond,h,ma[i])
        #print(i,pbond,hess_num[i][i],post)
        for j in range(i+1,post.size):  # off diagonals
            # uses d^2f(x)/ dx dy = (f(x+h,y+h)+f(x-h,y-h)-f(x+h,y-h)-f(x-h,y+h))/(4hh)
            tposi = post[i]
            tposj = post[j]
            post[i] = tposi - h
            post[j] = tposj - h
            pbondm1m1 = bond_force(bond_style,nbonds,bonds,bondcoeff,post.reshape(-1,3),acc_d)
            post[j] = tposj + h
            pbondm11 = bond_force(bond_style,nbonds,bonds,bondcoeff,post.reshape(-1,3),acc_d)
            post[i] = tposi + h
            pbond11 = bond_force(bond_style,nbonds,bonds,bondcoeff,post.reshape(-1,3),acc_d)
            post[j] = tposj - h
            pbond1m1 = bond_force(bond_style,nbonds,bonds,bondcoeff,post.reshape(-1,3),acc_d)
            post[i] = tposi
            post[j] = tposj

            hess_num[i][j] = (pbond11+pbondm1m1-pbondm11-pbond1m1)/(h*h*4.0)
            hess_num[i][j] /= math.sqrt(ma[i]*ma[j])
            hess_num[j][i] = hess_num[i][j]

#-------------------------------------------------
def bond(bond_style,nbonds,bonds,bondcoeff,pos,acc,masses):

    #print("Calculating bond forces")
    pbond = 0
    pbond = bond_force(bond_style,nbonds,bonds,bondcoeff,pos,acc)

    return pbond
    # check routines...
    #check_forces(bond_style,nbonds,bonds,bondcoeff,pos,acc,masses)
    #print(acc)
    #check_inm(bond_style,nbonds,bonds,bondcoeff,pos,acc,masses)
    #exit(1)

#----------------------------------------------------------
def bond_hess(bond_style,nbonds,bonds,bondcoeff,pos,masses,hessian):
    inm(bond_style,nbonds,bonds,bondcoeff,pos,masses,hessian)

#----------------------------------------------------------
def check_forces(bond_style,nbonds,bonds,bondcoeff,pos,acc,masses):
    #check forces

    acc.fill(0) # rezero forces
    pbond = bond_force(bond_style,nbonds,bonds,bondcoeff,pos,acc)

    acc_num = numpy.copy(acc)
    bond_num(bond_style,nbonds,bonds,bondcoeff,pos,acc_num)

    tol = 1e-6
    diff = acc_num-acc
    mv = max(diff.max(),abs(diff.min()))

    if(mv < tol):
        print("Forces Match, pbond =",pbond,mv)
    else:
        print("Forces DO NOT Match!")
        print("Analytical")
        print(acc)
        print("Numerical")
        print(acc_num)
        print("Diff = ",diff)

#----------------------------------------------------------
def check_inm(bond_style,nbonds,bonds,bondcoeff,pos,acc,masses):

    print("Calculating Hessian")
    hessian = numpy.zeros((pos.size,pos.size))
    bond_hess(bond_style,nbonds,bonds,bondcoeff,pos,masses,hessian)
    #print(hessian)

    #print(pos.size)
    hess_num = numpy.zeros((pos.size,pos.size))
    bond_hess_num(bond_style,nbonds,bonds,bondcoeff,pos,acc,masses,hess_num)
    #print(hess_num)

    hdiff = hess_num-hessian
    mv = max(hdiff.max(),abs(hdiff.min()))/hessian.max()

    rdiff = numpy.sum(numpy.sum(hdiff*hdiff,axis=1))
    tol = 1e-5
    print(rdiff,mv)
    if(rdiff<tol and mv < tol):
        pbond = bond_force(bond_style,nbonds,bonds,bondcoeff,pos,acc)
        print("Hessians Match, pbond =",pbond,mv,rdiff)
    else:
        print("Hessians DO NOT Match!")
        print(hessian)
        print(hess_num)
        print("Differences!")
        print("Diff = ",mv,rdiff,hdiff)
        print("Bummer!")
        exit(1)

    print("omega-squared")
    for i in range(nbonds):  # loop over bonds,
        itype = bonds[i][0]  # bond type
        mi = masses[bonds[i][1]][0]
        mj = masses[bonds[i][2]][0]
        mu =mi*mj/(mi+mj) # reduced mass
        if(bond_style==0): #Harmonic
            k0 = bondcoeff[itype][0] # use type to bond params
            r0 = bondcoeff[itype][1]
            print("Harmonic",2*k0/mu)
        elif(bond_style==1): #Morse:
            D = bondcoeff[bonds[i][0]][0] #idx D alpha r0
            alpha = bondcoeff[itype][1]
            r0 = bondcoeff[itype][2]
            ipos = pos[bonds[i][1]]
            jpos = pos[bonds[i][2]]
            dpos = jpos-ipos
            r =  math.sqrt(numpy.dot(dpos,dpos))
            dr = r-r0
            expar = math.exp(-alpha*dr)
            print("Morse",(2*alpha*alpha*D*(-1*expar+2*expar*expar))/mu)
        else: #??
            print("ERROR in bond_style")
            exit(1)

    w,v = numpy.linalg.eig(hessian)
    print("eigenvalus:",w)
    #print(v)

    w,v = numpy.linalg.eig(hess_num)
    print("eigenvalues num:", w)
    #print(v)

    print("Eigenvectors")
    for i in range(pos.size):
        print("")
        print(i,w[i],v[:,i])

    f = open("check.vmd","w")
    f.write("mol new\n")
    f.write("draw color blue\n")
    f.write("draw sphere {%f %f %f} radius .2\n" % (pos[0][0], pos[0][1], pos[0][2]))
    f.write("draw color red\n")
    f.write("draw sphere {%f %f %f} radius .2\n" % (pos[1][0], pos[1][1], pos[1][2]))
    for i in range(6):
        f.write("draw color %d\n" %(i+1))
        f.write("draw line {%f %f %f} {%f %f %f} width 3\n" % (pos[0][0], pos[0][1], pos[0][2], pos[0][0]+v[0][i],pos[0][1]+v[1][i], pos[0][2]+v[2][i]))
        f.write("draw line {%f %f %f} {%f %f %f} width 3\n" % (pos[1][0], pos[1][1], pos[1][2], pos[1][0]+v[3][i],pos[1][1]+v[4][i], pos[1][2]+v[5][i]))
    f.close()

    exit(1)
#endcheck

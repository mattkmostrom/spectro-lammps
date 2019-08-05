# output for mdpython

import numpy

#------------------------------------------------------------------
def write_dump_vel(velfile,istep,natoms,vel):
    global fowrite_vel

    try: fowrite_vel  # has fowrite been assigned yet?
    except NameError:
        fowrite_vel = open(velfile,"w")

    fowrite_vel.write("# {} step {}\n".format(natoms,istep))
    for i in range(natoms):
        line = str(vel[i][0]) + " " + str(vel[i][1]) + " " + str(vel[i][2]) + "\n"
        fowrite_vel.write(line)
    # fo.close()
    return 0

#------------------------------------------------------------------
def write_dump(dumpfile,istep,natoms,pos,aatype):
    global fowrite
    try: fowrite  # has fowrite been assigned yet?
    except NameError:
        fowrite = open(dumpfile,"w")

    fowrite.write("%d\n"% natoms)
    fowrite.write("# step %d\n" % istep)
    for i in range(natoms):
        if(aatype[i]==1):
            line = "O " + str(pos[i][0]) + " " + str(pos[i][1]) + " " + str(pos[i][2]) + "\n"
        elif (aatype[i]==2):
            line = "C " + str(pos[i][0]) + " " + str(pos[i][1]) + " " + str(pos[i][2]) + "\n"
        else :
            print("Error type not found?!")
            exit(1)
        fowrite.write(line)
    # fo.close()
    return 0

#---------------------------------------------------------------------------
def write_thermo(logfile,istep,natoms,masses,pos,vel,pot):

    global fothermo

    ke = .5*numpy.sum(masses*vel*vel)
    tpot = numpy.sum(pot)
    temp = 2.0*ke/(3.0*natoms)

    if(logfile==None):
        if(istep==0):
            print('#step, temp, etotal, ke, tpot, evdw, ebond')
#        print(istep,temp,ke+tpot,ke,tpot,pot[0],pot[1])
        print('{:8} {:8.5} {:8.5} {:8.5} {:8.5} {:8.5} {:8.5}'.format(istep,temp,(ke+tpot),ke,tpot,pot[0],pot[1]))
#        print('{:8} {:8.5} {:8.5} {:8.5} {:8.5} {:8.5} {:8.5}'.format(istep,temp,(ke+tpot)/natoms,ke/natoms,tpot/natoms,pot[0]/natoms,pot[1]/natoms))
    else:
        try: fothermo  # has fothermo been assigned yet?
        except NameError:
            fothermo = open(logfile,"w")
            fothermo.write('# step, temp, etotal, ke, tpot, evdw, ebond\n')

        line = str(istep)+" "+str(temp)+" "+str((ke+tpot))+" "+str(ke)+" "+str(tpot)+" "+str(pot[0])+" "+str(pot[1])+"\n"
        fothermo.write(line)

    # fo.close()
    return ke+tpot

#---------------------------------------------------------------------------
def write_init(initfile,istep,natoms,atypes,nbonds,tbonds,box,mass,pos,vel,bonds,aatype):

    foinit = open(initfile,"w")

    foinit.write("Lammps init file written by MDPython after %d\n"% istep)
    foinit.write("\n")
    foinit.write("%d atoms\n"% natoms)
    foinit.write("%d atom types\n"% atypes)
    foinit.write("%d bonds\n"% nbonds)
    foinit.write("%d bond types\n"% tbonds)
    foinit.write("\n")
    line = str(-box[0]/2) + " " + str(box[0]/2) + " xlo xhi\n"
    foinit.write(line)
    line = str(-box[1]/2) + " " + str(box[1]/2) + " ylo yhi\n"
    foinit.write(line)
    line = str(-box[2]/2) + " " + str(box[2]/2) + " zlo zhi\n"
    foinit.write(line)
    foinit.write("\n")
    foinit.write("Masses\n")
    foinit.write("\n")
    for i in range(atypes):
        line = str(i+1) + " " + str(mass[i]) + "\n"
        foinit.write(line)

    foinit.write("\n")
    foinit.write("Atoms # full style: index type group charge x y z")
    foinit.write("\n")
    for i in range(natoms):
        line = str(i+1) + " " + str(aatype[i]) + " 1 0 " + str(pos[i][0]) + " " + str(pos[i][1]) + " " + str(pos[i][2]) + "\n"
        foinit.write(line)

    foinit.write("\n")
    foinit.write("Velocities # vx vy vz\n")
    foinit.write("\n")
    for i in range(natoms):
        line = str(i+1) + " " + str(vel[i][0]) + " " + str(vel[i][1]) + " " + str(vel[i][2]) + "\n"
        foinit.write(line)

    foinit.write("\n")
    foinit.write("Bonds # type i j\n")
    foinit.write("\n")
    for i in range(nbonds):
        line = str(i+1) + " " + str(bonds[i][0]+1) + " " + str(bonds[i][1]+1) + " " + str(bonds[i][2]+1) + "\n"
        foinit.write(line)

    foinit.write("\n")
    foinit.close()
    return 0

def write_inm(istep,hessian):
    #print(istep,hessian)
    print(istep)
    # print("eigenvalus:",w)

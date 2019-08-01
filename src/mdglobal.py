

def globalinit():
    print("initializing global variables")

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
    global abtype       # array of bond types    
    global logfile      # file to output thermodata
    global hessian      # hessian matrix


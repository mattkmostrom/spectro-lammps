# import routines for MDPython
#TODO
#Fix replicate method
#Add angles

import re
import numpy as np

global natoms
global atypes
global nbonds
global tbonds
global box

box = np.zeros(3)

re_dict_data = {
        'natoms': re.compile(r'^(?P<natoms>\d+)\s+atoms\n'),
        'atypes': re.compile(r'^(?P<atypes>\d+)\s+atom types\n'),
        'nbonds': re.compile(r'^(?P<nbonds>\d+)\s+bonds\n'),
        'tbonds': re.compile(r'^(?P<tbonds>\d+)\s+bond types\n'),
        'box_x': re.compile(r'^(?P<box_x>[.\d-]+\s+[.\d-]+)\s+xlo\s+xhi\n'),
        'box_y': re.compile(r'^(?P<box_y>[.\d-]+\s+[.\d-]+)\s+ylo\s+yhi\n'),
        'box_z': re.compile(r'^(?P<box_z>[.\d-]+\s+[.\d-]+)\s+zlo\s+zhi\n')
        }

data_lst = [key for key in re_dict_data]

def parse_line(line,d):
    for key, rx in d.items():
        match = rx.search(line)
        if match:
            return key, match

    return None, None

def readinvals(datafile):
    with open(datafile,'r') as f:
        data = [-1]*7
        line = f.readline()
        while line:
            key, match = parse_line(line,re_dict_data)
            if match:
                val = match.group(key)
            if key == 'natoms':
                data[0] = int(val)
            if key == 'atypes':
                data[1] = int(val)
            if key == 'nbonds':
                data[2] = int(val)
            if key == 'tbonds':
                data[3] = int(val)
            if key == 'box_x':
                val = val.split()
                temp = float(val[1]) - float(val[0])
                data[4] = temp
            if key == 'box_y':
                val = val.split()
                temp = float(val[1]) - float(val[0])
                data[5] = temp
            if key == 'box_z':
                val = val.split()
                temp = float(val[1]) - float(val[0])
                data[6] = temp
            line = f.readline()
    
    for i in range(len(data)):
        if data[i] == -1:
            print('Not enough data found. Check init file for proper formatting.')
            print('Not found:', data_lst[i])
            exit(1)

    return data

re_dict_arrays = {
        'masses': re.compile(r'^Masses'),
        'atoms': re.compile(r'^Atoms'),
        'vels': re.compile(r'^Velocities'),
        'bonds': re.compile(r'^Bonds')
        }

def make_arrays(datafile,reps):
    global natoms, atypes, nbonds, tbonds, box
    
    natoms, atypes, nbonds, tbonds, box[0], box[1], box[2] = readinvals(datafile) 
    
    with open(datafile,'r') as f:
        mass = []
        aatype = []
        pos = []
        vel = []
        masses = []
        bonds = []
        
        line = f.readline()
        while line:
            key, match = parse_line(line,re_dict_arrays)
            count = 1
            if key == 'masses':
                line = f.readline()
                while len(line.split()) == 0:
                    line = f.readline()
                for i in range(atypes):
                    words = line.split()
                    if int(words[0]) != count:
                        print('Error while assigning masses')
                        print('Expecting',count,'got',words[0])
                        exit(1)
                    m = float(words[1])
                    mass.append(m)
                    count += 1
                    line = f.readline()
                print('Assigned masses')
            if key == 'atoms':
                line = f.readline()
                while len(line.split()) == 0:
                    line = f.readline()
                for i in range(natoms):
                    words = line.split()
                    if int(words[0]) != count:
                        print('Error while assigning atoms')
                        print('Expecting',count,'got',line.split()[0])
                        exit(1)
                    atype = int(words[1])
                    aatype.append(atype)
                    pos.append([float(words[j + 4]) for j in range(3)])
                    masses.append([mass[atype - 1] for j in range(3)])
                    count += 1
                    line = f.readline()
                print('Assigned atom types and positions')
            if key == 'vels':
                line = f.readline()
                while len(line.split()) == 0:
                    line = f.readline()
                for i in range(natoms):
                    words = line.split()
                    if int(words[0]) != count:
                        print('Error while assigning velocities')
                        print('Expecting',count,'got',line.split()[0])
                        exit(1)
                    vel.append([float(words[j + 1]) for j in range(3)])
                    count += 1
                    line = f.readline()
                print('Assigned velocities')
            if key == 'bonds':
                line = f.readline()
                while len(line.split()) == 0:
                    line = f.readline()
                for i in range(nbonds):
                    words = line.split()
                    if int(words[0]) != count:
                        print('Error while assigning velocities')
                        print('Expecting',count,'got',line.split()[0])
                        exit(1)
                    bonds.append([int(words[j + 1])-1 for j in range(3)])
                    count += 1
                    line = f.readline()
                print('Assigned bonds')
            line = f.readline()
    
    if len(vel) == 0:
        vel = np.zeros((natoms,3))
        print('Velocities will start at 0')
    if len(bonds) == 0:
        print('No bonds found')

    pos_copy = pos[:]
    vel_copy = vel[:]
    bonds_copy = bonds[:]
    count = 1 

    for i in range(reps[0]):
        for j in range(reps[1]):
            for k in range(reps[2]):
                if i == j == k == 0:
                    continue
                offset = [i*box[0],j*box[1],k*box[2]]
                for vec in pos_copy:
                    temp = [sum(x) for x in zip(vec,offset)]
                    pos.append(temp)
                for vec in vel_copy:
                    np.append(vel,vec)
                bond_offset = [0,count*len(pos_copy),count*len(pos_copy)]
                for bond in bonds_copy:
                    temp = [sum(x) for x in zip(bond,bond_offset)]
                    bonds.append(temp)
                count += 1

    fact = np.prod(reps)
    natoms *= fact
    nbonds *= fact
    box = [np.prod(x) for x in zip(box,reps)]

    return np.array(mass), np.array(aatype), np.array(pos), np.array(vel), np.array(masses), np.array(bonds)

re_dict_sysvals = {
        'nsteps': re.compile(r'^run\s+(?P<nsteps>\d+)'),
        'dt': re.compile(r'^timestep\s+(?P<dt>[\d.]+)'),
        'initfile': re.compile(r'^read_data\s+(?P<initfile>[_A-z.\d]+)'),
        'ithermo': re.compile(r'^thermo\s+(?P<ithermo>\d+)'),
        'dump': re.compile(r'^dump\s+[A-z]+\s+[A-z]+\s+[A-z]+\s+(?P<dump>\d+\s+[_A-z.\d]+)'),
        'bond_style': re.compile(r'^bond_style\s+(?P<bond_style>[A-z]+)'),
        'logfile': re.compile(r'^log\s+(?P<logfile>[_A-z.\d]+)'),
        'inm': re.compile(r'^#inm\s+(?P<inm>[_A-z.\d]+\s+\d+)'),
        'reps': re.compile(r'^replicate\s+(?P<reps>\d+\s+\d+\s+\d+)'),
        'fix': re.compile(r'^fix\s+[_A-z\d]+\s+[_A-z\d]+\s+(?P<fix>[A-z]+)')
        }

sysvals_lst = ['dt','initfile','bond_style','idump','dumpfile','ithermo','logfile','inmfile','inmo','nsteps']

def readsysvals(infile):
    with open(infile,'r') as f:
        data = [-1]*10
        bondcoeff = []
        reps = [1,1,1]
        line = f.readline()
        while line:
            key, match = parse_line(line,re_dict_sysvals)
            if match:
                val = match.group(key)
            if key == 'nsteps':
                data[9] = int(val)
            if key == 'dt':
                data[0] = float(val)
            if key == 'initfile':
                initfile = val.strip(' ')
                data[1] = initfile
                natoms, atypes, nbonds, tbonds, box[0], box[1], box[2] = readinvals(initfile) 
            if key == 'ithermo':
                data[5] = int(val)
            if key == 'logfile':
                data[6] = val.strip(' ')
            if key == 'dump':
                val = val.split()
                data[3] = int(val[0])
                data[4] = val[1]
            if key == 'inm': 
                val = match.group(key).split()
                data[7] = val[0]
                data[8] = int(val[1])
            if key == 'bond_style':
                data[2] = val.strip(' ')
                if re.search('harmonic',val,flags=re.IGNORECASE):
                    bond_style = 0
                    print('Reading in harmonic bond coefficients for',tbonds,'types')
                    for i in range(tbonds):
                        line = f.readline()
                        bondcoeff.append([float(line.split()[j+2]) for j in range(2)])
                elif re.search('morse',val,flags=re.IGNORECASE):
                    bond_style = 1
                    print('Reading in morse bond coefficients for',tbonds,'types')
                    for i in range(tbonds):
                        line = f.readline()
                        bondcoeff.append([float(line.split()[j+2]) for j in range(3)])
                else:
                    print('Unrecognized bond type')
            if key == 'fix':
                fix_type = val
                words = line.split()
                var_lst = words[4:]
            if key == 'reps':
                vals = match.group(key).split()
                reps = [int(num) for num in vals]
            line = f.readline()
    
    bondcoeff = np.array(bondcoeff)

    for i in range(len(data)):
        if data[i] == -1:
            print('Not enough data found. Check init file for proper formatting.')
            print('Not found:', sysvals_lst[i])
            exit(1)

    return data, bond_style, bondcoeff, reps, fix_type, var_lst

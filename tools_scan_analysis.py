#!/usr/bin/python
#nada
import numpy as np
from numpy.linalg import *
import os
import glob

#def parseALL(filename,first,second):
def parseALL(filename,first,second):
    parsing = False
    name = []
    f = open(filename)
    for line in f.readlines():
        if line.startswith(first):
            parsing = True
            continue
        elif line.startswith(second):
            parsing = False
        if parsing:
            a = line.split()
            name.extend(a)
#    name = map(float,name)
    return name

#################################################
#################################################
#################################################

def parser(filename,first,second):
    whole_data = []
    parse = False
    with open(filename,'r') as out:
     for line in out:
        if str(first) in line:
            parse = True
#            continue
        elif str(second) in line:
            parse = False
        if parse:
            whole_data.append(line)
    out.close()
    return whole_data

##############################
## Parser for nested list:
##############################

def parser_lists(data,first,second):
    whole_data = []
    parse = False
    for line in data:
        if str(first) in line:
            parse = True
        elif str(second) in line:
            parse = False
        if parse:
            whole_data.append(line)
    return whole_data



def chunks_optim(filename):
    with open(filename) as f:
        lines = f.read().splitlines()
    input_or = [i for i,k in enumerate(lines) if 'Input orientation:' in k]
    chunks = []
    a = 0
    for ind in input_or:
        chunks.append(lines[a:ind])
        a = ind
        if input_or.index(ind) == len(input_or)-1:
            chunks.append(lines[ind:])
        
    b = []
    for geom in chunks:
        for s in geom:
            if 'Optimization completed' in s:
             b.append(geom)

    return (b)

def extract_distance(filename,atom1,atom2):
    data = chunks_optim(filename)
    all_dist = []
    for i in data:
        a = parser_lists(i,'Input orientation:','Distance matrix')[5:-1]
        b = [ i.split()[3:] for i in a]
        c = np.array([map(float,i) for i in b])
        dist = norm(c[atom2-1] - c[atom1-1])
        all_dist.append(dist)
    all_F1, all_F2 = [],[]
    for j in data:
        ad = parser_lists(j,'(Enter /home/apps/Chem/g09/l716.exe)','Cartesian Forces')[5:-1]
        bd = [ i.split()[2:] for i in ad]
        cd = np.array([map(float,i) for i in bd])
        to_nN = 82.387
        F1 = norm(cd[atom1-1]) *to_nN
        F2 = norm(cd[atom2-1]) *to_nN
        all_F1.append(F1)
        all_F2.append(F2)
    return (all_dist, all_F1, all_F2 )

def extract_E(filename):
    data = chunks_optim(filename)
    all_E = []
    for i in data:
        a = parser_lists(i,')     EIGENVALUE','Final one electron')
        b = [ k.split() for k in a if 'EIGENVALUE' in k]
        c = map(float,[j[-1] for j in b])
        all_E.append(c)
    return all_E
    
# Test type of calculation:
def extract_E_W(filename,atom1,atom2):
    dist0 = extract_distance(filename,atom1,atom2)[0][0]
    dists = extract_distance(filename,atom1,atom2)[0]
    dif_dists = np.array([i-dist0 for i in dists]) * 1.889725989
    forces = np.array(extract_distance(filename,atom1,atom2)[1]) /82.387
    W = [a*b for a,b in zip(forces,dif_dists)]
    E = extract_E(filename)
    E_W = [a+np.array(b) for a,b in zip(W,E)]
    E_W1 = [ a.tolist() for a in E_W]
    return E_W1
    

def writer_dist_F_E_CASSCF(filename,atom1,atom2,fileplot):
    dist_list = extract_distance(filename,atom1,atom2)[0]
    F1_list = extract_distance(filename,atom1,atom2)[1]
    E_list = extract_E(filename)
    EW_list = extract_E_W(filename,atom1,atom2)
    all_dat = zip(dist_list,F1_list,E_list,EW_list)
    with open(fileplot,'w') as f:
        for i in all_dat:
           a = map(str,i[2])
           b = map(str,i[3])
           f.write('%s  %s  %s  %s \n' % (i[0], i[1], '  '.join(a), '  '.join(b) ))
    return 



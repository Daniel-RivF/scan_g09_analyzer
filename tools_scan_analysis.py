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
    

def writer_dist_F_E_CASSCF(filename,atom1,atom2,fileplot):
    dist_list = extract_distance(filename,atom1,atom2)[0]
    F1_list = extract_distance(filename,atom1,atom2)[1]
    E_list = extract_E(filename)
    all_dat = zip(dist_list,F1_list,E_list)
    with open(fileplot,'w') as f:
        for i in all_dat:
           a = map(str,i[2])
           f.write('%s  %s  %s  \n' % (i[0], i[1], '  '.join(a) ))
    return 


###################################################################################
###################################################################################
# This function extracts times ana an  user-defined distance from a molcas MD-Tully output.
###################################################################################
###################################################################################

def getdistances2(filename,Natom,atom1,atom2):
    whole_data = parser(filename,' New Coordinates (','--- Stop Module')
    #xyz in bohr
    tofs = 2.418884326505e-2
    bohrtoA = 0.529177249
    aa = [ j for j in whole_data if "---------------" not in j ]
    ab = [ k for k in aa if "New Coordinates " not in k ]
    ac = [l.split()  for l in ab]
    ad = [i for i in ac if i!=[] ]
    ae = [i[2:5] for i in ad if 'X' not in i]
    af = [map(float,m) for m in ae]
    ag =  [ [float(y)*bohrtoA for y in x] for x in af]
    xyz = [ag[h:h+Natom] for h in range(0, len(ag), Natom)]
    distances = [norm(at[atom1-1]-at[atom2-1]) for at in np.array(xyz)]
    #times
    t =  parser(filename,'New Coordinates (time=','------') 
    tt = [l.split()  for l in t]
    ttt = [i[3] for i in tt] 
    tim = [s.strip('a.u.):') for s in ttt]
    times = tofs * np.array(map(float,tim))
    # HOPPING INFO
    time_hop = parser(filename,'A HOP event is detected','No.')
    a = [ j for j in time_hop if "time=" in j]
    b = [l.split()  for l in a]
    c = [a[3::] for a in b]
    cc = [i for sublist in c for i in sublist]
    d  = [ ''.join(i.split())[:-6] for i in cc]
    times_hop = np.array(map(float,d)) * tofs
    ah  = [tim.index(i) for i in d]
    disthop = [distances[i] for i in ah]
    return (distances[1::],times.tolist()[1::],disthop,times_hop)


#################################################################################
# This function extracts the Energies and populations from Tully-Molcas MD output.
# The hop point and time too
##################################################################################
# Ojo con listt
# listt,listd,times_hop,disthop from previous function:
def projectVelocityHOP(filename,atom1,atom2):
    atoms=parser(filename,'Coordinates and Masses of Atoms, in au and A',' The Moment of Inertia Tensor ')
    Natom = len(atoms) - 3
    listv = parser(filename,'Velocities before Hop:','Velocities after Hop:')
    listv1 = [ i.split() for i in  listv[1:Natom + 1] ]
    listv2 = [map(float,i) for i in listv1]
    if listv != []:
       vatom1 = np.array(listv2[atom1-1])
       vatom2 = np.array(listv2[atom2-1])
       # Find index of first of first hop
       op = parser(filename,'A HOP event is','Molecular Dynamics specifications')
       indice = [i for i, x in enumerate(op) if 'Old Coordinates' in x][0]
       geom = op[indice+4:indice+Natom+4]
       geom1 =  [i.split() for i in geom]
       geom2 = [map(float,i[2:5]) for i in geom1]
       xatom1 = np.array(geom2[atom1-1])
       xatom2 = np.array(geom2[atom2-1])
       line21 = (xatom2 - xatom1)/norm(xatom2 - xatom1)
       proj1_21 = np.dot(vatom1,line21)
       proj2_21 = np.dot(vatom2,line21)
       direction1 = 0
       direction2 = 0
       if proj1_21 < 0 :
          direction1 = 1
       else :
           direction1 = -1
       if proj2_21 < 0 : 
           direction2 = -1
       else:
           direction2 = 1
    # direction_i = 1  ==> contraction (projection has same direction that the defined 21 line)
    # direction_i = -1 ==> extension (projection has opposite direction that the defined 21 line)
       return (direction1,direction2)
    else:
       return ('no hop')
    

#    listv2 = [item for sublist in listv1 for item in sublist]
#    listv3 = map(float,listv2)

    



def getEnersPop(filename,listt):
#    times = getdistances2(filename,Natom,atom1,atom2)[1]
#    disthop = getdistances(filename,Natom,atom1,atom2)[3]
    whole = parser(filename,'OOL','------')
    a = [ i.split() for i in whole]
    # Filter when HOP ALLOWED
    hop_E = []
    checkhop = [ 'ALLOWED' in n for n in a]
    indexes = [ i for i, x in enumerate(checkhop) if x == True ]
    if indexes != [] :
        for i in indexes:
            hop_point = a[i-1]
            hop_E.append(hop_point)
    else:
        hop_E = 'No hop in %s' % filename
   #Extract indexes with no data OOLGnuplot:
    nodata = ['OOLgnuplt:' not in i for i in a] 
    n = [ i for i, x in enumerate(nodata) if x == True ]
    if n != []:
       filtered = [l for k, l in enumerate(a) if k not in n]   
       aa = [ m[1::] for m in filtered ]
       datas = [ map(float,i) for i in aa]
    else:
        aa = [ m[1::] for m in a]
        datas = [ map(float,i) for i in aa]
        times_hop = 'no hop time'
    if len(datas) == len(listt) - 1:
       del listt[-1]
    return (datas, hop_E, listt)

######## VAYA PEDAZO DE MIERDA ############


def joinfiles(filename1,filename2):
    a = []
    with  open(filename1,'r') as f1:
     for line in f1:
       a.append(line)
    aa = [l.split()  for l in a]
    l1 =  [ map(float,i) for i in aa]
    last_t = map(float,aa[-1])[0]
    b = [] 
    with  open(filename2,'r') as f2:
      for line in f2:
        b.append(line)
    bb = [l.split()  for l in b]
    l2 = [ map(float,i) for i in bb] 
    final = [ [i[0] + last_t] + i[1:] for i in l2]
    total = l1 + final
   # f = open()
#    fileout = os.path.splitext(filename1)[0]
#    f = open(fileout + 'F', 'w')
#    n = [map(str,n) for n in l1+final]
#    # Writing a nested list to a file:
#    f.writelines(' '.join(i) + '\n' for i in n)
#    f.close()
    return total


########################################################
################### Print all info here #################
########################################################


def writer_dat(listE,listt,listd,fileout='EnersPoP.txt'):
    if [len(listE)]*2 == [len(listt)-1,len(listd)-1]:
       del listt[-1]
       del listd[-1]
    f=open(fileout,'w')
    for i in zip(listE,listt,listd):
        a = map(str,i[0])
        f.write('%s %s %s' % (i[1], i[2], '  '.join(a) ))
        f.write("\n")
    return

################################################
# Extracted lists from main ( pointhop, timehop).
################################################

def writerhop(fileout, pointhop, timehop, disthop):
    f = open(fileout,'w')
    for i in zip(pointhop,timehop,disthop):
        f.write('%s %s %s' % (i[1], i[2], '  '.join(i[0][1::]) ))
        f.write("\n")
    return 








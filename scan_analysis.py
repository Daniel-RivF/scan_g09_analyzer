#!/usr/bin/python
import tools_scan_analysis as tools
filename = raw_input("Name of the g09 scan output file:   ")
atom1 = int(raw_input("1st pulling point (atom label):    "))
atom2 = int(raw_input("2nd pulling point (atom label):    "))
type_job = int(raw_input("Type of calculation  \n (Type 1 for CASSCF and 2 for DFT):    " ))
fileplot = raw_input(" name of the output file   :" )

if type_job == 1:
    tools.writer_dist_F_E_CASSCF(filename,atom1,atom2,fileplot)
#elif type_job == 2:
#    tools.writer_dist_F_E_DFT(filename,atom1,atom2,fileplot)


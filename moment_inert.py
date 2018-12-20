#! /usr/bin/python

import sys, os
from optparse import OptionParser
from math import exp, log, expm1, log1p, pi, pow, sqrt
import numpy as np
import matplotlib.pyplot as plt

homedir = "/home/stvogt/"

h = 6.62606957e-34
R = 8.3144621
Rkcal = 1.987191683e-3
T = 298.15
#T = 100.
k = 1.3806488e-23
c = 29979245800
alpha = 4
w0 = 100
mass = {"C" : 12.0107, "H": 1.00794, "O":15.9994, "N": 14.00674, "Al":26.981538}
bohr2angs = 0.52917721092


def read_coordsG09(output):
    atom_list = []
    x_coord = []
    y_coord = []
    z_coord = []
    try:
        lines = open(output,'r').readlines()
    except IOError:
        print "\n\nError!!   No output found in:  "+os.getcwd()
        sys.exit(1)

    for lineNum in range(0,len(lines)):
        line = lines[lineNum]
        if "Charge =" in line:
            for lineNum1 in range(lineNum+1,len(lines)):
                line1 = lines[lineNum1]
                if len(line1.split()) > 3:
                    atom_list.append(line1.split()[0])
                    x_coord.append(float(line1.split()[1]))
                    y_coord.append(float(line1.split()[2]))
                    z_coord.append(float(line1.split()[3]))
                else:
                    break
    #print "Input Coordinates: "
    #for i in range(0,len(atom_list)):
    #    print atom_list[i] + "    "+str(x_coord[i])+"     "+str(y_coord[i])+"     "+str(z_coord[i])
    #print '\n'
    return (atom_list, x_coord, y_coord, z_coord)

def read_coords_orca(output):
    atom_list = []
    x_coord = []
    y_coord = []
    z_coord = []
    try:
        lines = open(output,'r').readlines()
    except IOError:
        print "\n\nError!!   No output found in:  "+os.getcwd()
        sys.exit(1)

    for lineNum in range(0,len(lines)):
        line = lines[lineNum]
        if "CARTESIAN COORDINATES (ANGSTROEM)" in line:
            for lineNum1 in range(lineNum+2,len(lines)):
                line1 = lines[lineNum1]
                if len(line1.split()) == 4:
                    atom_list.append(line1.split()[0])
                    x_coord.append(float(line1.split()[1]))
                    y_coord.append(float(line1.split()[2]))
                    z_coord.append(float(line1.split()[3]))
                else:
                    break
    #print "Input Coordinates: "
    #for i in range(0,len(atom_list)):
    #    print atom_list[i] + "    "+str(x_coord[i])+"     "+str(y_coord[i])+"     "+str(z_coord[i])
    #print '\n'
    return (atom_list, x_coord, y_coord, z_coord)

def center_of_mass(coords):
    X_coord = []
    Y_coord = []
    Z_coord = []
    Xcm_tmp=0.0
    Ycm_tmp=0.0
    Zcm_tmp=0.0
    tot_mass = 0.0
    for i in range(0,len(coords[0])):
        tot_mass = tot_mass + mass[coords[0][i]]
    for i in range(0,len(coords[0])):
        Xcm_tmp = Xcm_tmp + mass[coords[0][i]]*coords[1][i]
        Ycm_tmp = Ycm_tmp + mass[coords[0][i]]*coords[2][i]
        Zcm_tmp = Zcm_tmp + mass[coords[0][i]]*coords[3][i]

    #print "Total Mass: "+str(tot_mass)+" amu\n"
    Xcm = Xcm_tmp/tot_mass  
    Ycm = Ycm_tmp/tot_mass
    Zcm = Zcm_tmp/tot_mass

    for i in range(0,len(coords[0])):
        X_coord.append(coords[1][i]-Xcm)
        Y_coord.append(coords[2][i]-Ycm)
        Z_coord.append(coords[3][i]-Zcm)
    #print "Center of Mass:  "+ str(Xcm) +'  '+ str(Ycm) + '  ' + str(Zcm)+'\n'

    #print "Center of Mass Coordinates: "
    #for i in range(0,len(coords[0])):
    #    print coords[0][i] + "    "+str(X_coord[i])+"     "+str(Y_coord[i])+"     "+str(Z_coord[i])
    #print '\n'
 
    return (coords[0],X_coord,Y_coord,Z_coord)


def moment_inertia(coords):
    Ixx = 0.0
    Iyy = 0.0
    Izz = 0.0
    Ixy = 0.0
    Ixz = 0.0
    Iyz = 0.0
    for i in range(0,len(coords[0])):
        #print mass[coords[0][i]]
        #Diagonal
        Ixx = Ixx + mass[coords[0][i]]*(pow(coords[2][i],2)+ pow(coords[3][i],2))
        Iyy = Iyy + mass[coords[0][i]]*(pow(coords[1][i],2)+ pow(coords[3][i],2))
        Izz = Izz + mass[coords[0][i]]*(pow(coords[1][i],2)+ pow(coords[2][i],2))
        #Off Diagonal
        Ixy = Ixy - mass[coords[0][i]]*coords[1][i]*coords[2][i]
        Ixz = Ixz - mass[coords[0][i]]*coords[1][i]*coords[3][i]
        Iyz = Iyz - mass[coords[0][i]]*coords[2][i]*coords[3][i]
    Itens = [[Ixx,Ixy,Ixz],[Ixy,Iyy,Iyz],[Ixz,Iyz,Izz]]
    eigvals, eigvecs = np.linalg.eig(Itens)
    #print eigvals
    #print eigvecs
    return eigvals


#Main

#output = sys.argv[1]
#coords = read_coords_orca(output)
##print coords
#COM_coords = center_of_mass(coords)
#
#moment_inertia(COM_coords)








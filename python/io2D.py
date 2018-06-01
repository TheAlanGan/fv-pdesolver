# -*- coding: utf-8 -*-
"""
Alan Gan - 11/29/17
M578 - 2D Parallel

Inputs/Outputs
"""
import numpy as np

def getInput(filename):
    filename = filename        # The name of file containing input values

    with open(filename) as f:        # Opens the txt file
        inputs = f.read().splitlines()    # Converts the txt file into array

    count = 0
    with open(filename) as f:
        for line in f:
            line = line.partition('#')[0]  # Strips everything after '#'
            line = line.rstrip()
            inputs[count] = line
            count += 1


    MMr = int(inputs[0])
    MMz = int(inputs[1])
    tEnd = float(inputs[2])
    dtout = float(inputs[3])
    factor = float(inputs[4])
    RIn = float(inputs[5])
    ROut = float(inputs[6])
    ZIn = 0.0
    ZOut = np.pi
    D = float(inputs[7])
    nodeCountType = bool(eval(inputs[8]))
    # nodeCountType (nct): True if input is MMr/MMz . False if input is Mr/Mz

    return MMr, MMz, tEnd, dtout, factor, RIn, ROut, ZIn, ZOut, D, nodeCountType
#=====================================================================

def otherValues(MMr, MMz, tEnd, dtout, factor, RIn, ROut, ZIn, ZOut, D, nct):
    dr = (ROut - RIn) / MMr
    dz = (ZOut - ZIn) / MMz
    Mr = MMr
    Mz = MMz

    if nct == True:
        dr = 1.0/MMr
        dz = 1.0/MMz
        Mr = int((ROut-RIn)*MMr)
        Mz = int((ZOut - ZIn)*MMz)

    dtEXPL = 1.0/( (2*D/dr**2) + (2*D/dz**2) )           # dt based on CFL condition
    dt = factor * dtEXPL           # Size of time-steps
    MaxSteps = int(tEnd/dt) + 1    # Number of time-steps

    return dr, dz, Mr, Mz, dtEXPL, dt, MaxSteps

#=====================================================================
def COMPARE(D, Mr, Mz, nsteps, time, output, R, Z, U): # For debugging

    # Creates the output file
    values = open('o.'+output, 'w')

    # Initiates the Uex Array
    Uex = np.zeros((Mz+2, Mr+2))

    # Initiating the Maximum Error Variable
    ERRmax = 0.0

    """_____Computing Exact Values and Comparing_____"""
    for i in range(0, Mr+2):
        for j in range (0, Mz+2):
            Uex[j, i] = np.exp(-time) * np.log(R[i]) * np.sin(Z[j])

            ERRij = abs(U[j, i] - Uex[j, i])#/Uex[j, i]#*100 #Remove other '#' for relative error in percent
            ERRmax  = max(ERRij, ERRmax)    # ERRmax is max error across all i's
    """_____END OF COMPARISON_____"""


    """_____THE OUTPUTTING TO TEXT FILE_____"""
    values.write("# U profile at t = " + str(time) + ";  nsteps = " + str(nsteps) +"\n")
    values.write("# Largest Error at this time: " + str(ERRmax) +"\n")
    values.write('{0:<13} {1:>18} {2:>18} {3:>18} {4:>18}'.format('# Z[j]', 'R[i]', 'U[i]', 'Uex[i]', 'Error') + "\n")

    for j in range(0, Mz+2):
        for i in range (0, Mr+2):
            values.write('{0:<13} {1:>18} {2:>18} {3:>18} {4:>18}'.format(Z[j], R[i], U[j,i], Uex[j,i], abs(U[j, i]-Uex[j, i]) )+ "\n")
    values.write("\n")
    """_____END OF OUTPUTTING_____"""

    # Closes the output file. It can now be read.
    values.close()

    # Returns the maximum error over all Control Volumes
    return ERRmax

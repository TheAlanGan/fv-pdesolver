# -*- coding: utf-8 -*-
"""
Alan Gan - 11/29/17
M578 - 2D Parallel
main master

"""
import sys
import numpy as np
from mpi4py import MPI

import io2D as io
import setup2D as setup
import update2D as update
import messaging2D as ms

def MASTER(nWRs, myID):
    try:
        filename = str(sys.argv[1])
    except IndexError:
        filename = 'inputs.txt'
    ### Looks for name of input file in execution arguments
    ### If no file name listed, then it looks for 'inputs.txt'
    ###     in the same directory
    ### For Example:
    ###         $     mpiexec -n 3 python main2D.py "MM32.txt"

    try:
    	outputName = str(sys.argv[2])
    except IndexError:
    	outputName = 'values'
    ### This does the same thing as above. If nothing in execution arguments
    ### Then the default output name is  o.values##.txt  where ## is a number.

    outputExtension = '.txt'
    numOut = 0 # To differentiate each output. Helps with outputting
    output = outputName + str(numOut) + outputExtension
    numOut += 1


    """_____GETTING PARAMETERS_____"""
    # Gets the parameters from input (text) file
    MMr, MMz, tEnd, dtout, factor, RIn, ROut, ZIn, ZOut, D, nct = io.getInput(filename) # nct is NodeCountType (see io2D.py)

    # Calculates the other necessary parameters based on above values.
    dr, dz, Mr, Mz, dtEXPL, dt, MaxSteps = io.otherValues(MMr, MMz, tEnd, dtout, factor, RIn, ROut, ZIn, ZOut, D, nct)
    """_____END OF PARAMETERS_____"""


    """_____PACKING AND SENDING PARAMETERS TO WORKERS_____"""
    parms = np.empty(11, dtype=np.float64) # Float parameters go into 'parms'
    parms[0] = tEnd
    parms[1] = dtout
    parms[2] = factor
    parms[3] = RIn
    parms[4] = ROut
    parms[5] = ZIn
    parms[6] = ZOut
    parms[7] = D
    parms[8] = dr
    parms[9] = dz
    parms[10] = dt

    iparms = np.empty(3, dtype=np.int)  # Integer parameters go into 'iparms'
    iparms[0] = Mr
    iparms[1] = Mz
    iparms[2] = MaxSteps

    MPI.COMM_WORLD.Bcast([parms, MPI.DOUBLE], root=0) # Broadcasts the parameters to workers
    MPI.COMM_WORLD.Bcast([iparms, MPI.INT], root=0) # Broadcasts the int parameters to workers


    """_____END OF PACKING AND SENDING_____"""


    """_____THE INITIALIZATION_____"""
    time = 0.0
    tout = dtout
    ERRmax = 0.0

    # Gets the Mesh and area arrays (So the master can output)
    R, Z, Ar, Az = setup.MESH(RIn, ROut, ZIn, ZOut, dr, dz, Mr, Mz, myID, nWRs)        # Constructs Mesh

    # Sets up Flux Arrays and Initial Conditions in U array (So the master can output)
    U, Fr, Fz = setup.INIT(time, R, Z, Mr, Mz)        # Initiates the U and F arrays
    """_____END OF INITIALIZATION_____"""
    

    U = ms.RECV_output_MPI(nWRs, Mr, Mz, U)
    ERRmax = io.COMPARE(D, Mr, Mz, 0, 0, output, R, Z, U) # For debugging

#The Time-stepping
    for nsteps in range(1, MaxSteps + 1):

        MPI.COMM_WORLD.Barrier()

        """_____THE UPDATING_____"""
        time = nsteps * dt
        """_____END OF UPDATING_____"""


        """_____RECIEVING FROM WORKERS AND OUTPUTTING_____"""
        if (time >= dtout):# or ((time >= .1) and (time < .1+dt)):    #if time >= tout: #outputs according to Tout

            output = outputName + str(numOut) + outputExtension # The name of output file (each time is different)

            U = ms.RECV_output_MPI(nWRs, Mr, Mz, U)

            ERRmax = io.COMPARE(D, Mr, Mz, nsteps, time, output, R, Z, U) # For debugging
            tout = tout + dtout
            numOut += 1 # So the next output will output to a different file

            print "\ntime = ", time
            print "Max Error = ", ERRmax
        """_____END OF OUTPUTTING_____"""


        """_____EXECUTED AT END OF TIMESTEPPING_____"""
        if time >= tEnd:
            print "\nMMr = ", MMr, " MMz = ", MMz
            print "Mr = ", Mr, " Mz = ", Mz
            print "Done at time =", time, "\nafter nsteps =", nsteps
            print "Max Error =", ERRmax     # For Debugging
        """_____END OF THAT_____"""


#End of Time-stepping

    if time < tEnd:
        print "..... something went wrong ....."



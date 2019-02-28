# -*- coding: utf-8 -*-
"""
Alan Gan - 11/29/17
M578 - 2D Parallel
main worker

"""
import numpy as np
from mpi4py import MPI

import io2D as io
import setup2D as setup
import update2D as update
import messaging2D as ms

def WORKER(nWRs, myID):

    """_____RECIEVING PARAMETERS FROM MASTER_____"""
    parms = np.empty(11, dtype=np.float64)
    iparms = np.empty(3, dtype=np.int)

    MPI.COMM_WORLD.Bcast(parms, root=0)
    MPI.COMM_WORLD.Bcast(iparms, root=0)
    """_____END OF RECIEVING PARAMETERS_____"""



    """_____UNPACKING PARAMETER ARRAYS_____"""
    tEnd = parms[0]
    dtout = parms[1]
    factor = parms[2]
    RIn = parms[3]
    ROut = parms[4]
    ZIn = parms[5]
    ZOut = parms[6]
    D = parms[7]
    dr = parms[8]
    dz = parms[9]
    dt = parms[10]

    Mr = iparms[0]
    Mz = iparms[1]
    MaxSteps = iparms[2]
    """_____END OF UNPACKING_____"""


    """_____CONVERTING PARAMETERS INTO LOCAL PARAMETERS_____"""
    Mz = Mz / nWRs    # Only decomposing domain in Z direction

    # Dividing domain acording to which process
    ZInL = ZIn + dz * (myID - 1) * Mz
    ZOutL = ZInL + dz * Mz

    if (myID == nWRs):
        ZOutL = ZOut


    Me = myID
    NodeUP = Me + 1
    NodeDN = Me - 1
    """_____END OF CONVERSION_____"""



    """_____THE INITIALIZATION_____"""
    time = 0.0
    tout = dtout

    # Gets the Mesh and area arrays
    R, Z, Ar, Az = setup.MESH(RIn, ROut, ZInL, ZOutL, dr, dz, Mr, Mz, Me, nWRs)        # Constructs Mesh

    # Sets up Flux Arrays and Initial Conditions in U array
    U, Fr, Fz = setup.INIT(time, R, Z, Mr, Mz)        # Initiates the U and F arrays
    """_____END OF INITIALIZATION_____"""


    ms.SEND_output_MPI(Me, NodeUP, NodeDN, Mr, Mz, U)

#The Time-stepping
    for nsteps in range(1, MaxSteps + 1):

        MPI.COMM_WORLD.Barrier() # To keep all nodes/processes in sync

        U = ms.EXCHANGE_bry_MPI(nWRs, Me, NodeUP, NodeDN, Mr, Mz, U) # Exchanging boundaries with other nodes/processes

        """_____THE UPDATING_____"""
        # Need to update time first b/c Boundary Conditions are time dependant
        time = nsteps * dt
        U, Fr, Fz = update.FLUX(RIn, ROut, ZIn, ZOut, R, Z, D, dr, dz, Mr, Mz, Fr, Fz, U, time, Me, nWRs)
        U = update.PDE(ZIn, ZOut, R, Z, Ar, Az, dr, dz, Mr, Mz, dt, Fr, Fz, U, time, Me, nWRs)

        """_____END OF UPDATING_____"""


        """_____SENDING DATA TO MASTER_____"""
        if (time >= dtout):    #if time >= tout: #outputs according to Tout
            ms.SEND_output_MPI(Me, NodeUP, NodeDN, Mr, Mz, U)
            tout = tout + dtout
        """_____END OF SENDING_____"""


#End of Time-stepping

    if time < tEnd:
        print "..... something went wrong ..... myID = ", Me



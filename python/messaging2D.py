"""
Alan Gan - 11/29/17
M578 - 2D Parallel

2D parallel - messaging

This file will contain:
    RECV_output_MPI()
        - for the master to receive the values at time = dtout

    SEND_output_MPI()
        - for every worker to send the values at time = dtout

    EXCHANGE_bry_MPI()
        - for workers to share bdry values with neighbors
"""
from mpi4py import MPI
import numpy as np

#==========================================================================
def RECV_output_MPI(nWRs, Mr, Mz, U):

#    h = MPI.INT.Create_contiguous(x2)
#    MPI.INT.Commit()
#    h.Commit()   # might be this not entirely sure

#    x2 = M + 2
    localMz = Mz/nWRs
    collectiveU = np.zeros((Mz+2, Mr+2)) # The U array that contains the entire system

    count = 0  # For combining the U arrays from workers


    """ This next part receives the arrays from worker processess, strips the
        local boundaries (only keeping the REAL boundaries), and combines them
        into one array (the REAL mesh) which is sent to the master for output.
    """
    for i in range (1, nWRs + 1):
        buff = np.zeros((localMz + 2, Mr+2))  # buff is the memory buffer for the recieving array
#        Jme = (i-1)*M + 1    # Not sure what JME is for... ask prof
        msgtag = 1000 + i
        MPI.COMM_WORLD.Recv(buff, source = i, tag = msgtag )
#        print buff, "RECIEVED FROM ", i, '\n'

#        if i == 2:
#            print buff, "*****BEFORE DELETING"

        if (i != 1 and i!= nWRs): # If middle parts of mesh, deletes boundary values (top and bott)
            buff = np.delete(buff, [0,localMz+1], 0 )

        if (i == 1) and (i != nWRs): # If bottom part, deletes top boundary
            buff = np.delete(buff, [localMz+1], 0 )

        if (i == nWRs) and (i != 1): # If top part, deletes bottom boundary
            buff = np.delete(buff, [0], 0 )

        for x in range(0, buff[:,0].size) :
            collectiveU[count] = buff[x]
            count +=1
    """ END of the section referred to in the above comment.
    """

#    h.Free()
#    print collectiveU, "AFTER BEING PUT TOGETHER"
    return collectiveU

#===========================================================================

def SEND_output_MPI(Me, NodeUP, NodeDN, Mr, Mz, U):

    """ This method is run on worker processes and sends their U array to the
        master for output. Note: the master sorts out the boundaries.
    """

#    x2 = M + 2     # Size of local_U arrays

#    h = MPI.INT.Create_contiguous(x2) # Not sure what h is, but it works
#    MPI.INT.Commit()
#    h.Commit()    # Might be this, not sure.

    mster = 0
    msgtag = 1000 + Me

    MPI.COMM_WORLD.Send(U, dest = mster, tag = msgtag)

#    h.Free() # h is from the Create_contiguous... not sure why this runs

    return

#===========================================================================

def EXCHANGE_bry_MPI(nWRs, Me, NodeUP, NodeDN, Mr, Mz, U):

    h = MPI.DOUBLE.Create_contiguous(Mr+2)
    h.Commit()   # might be this not entirely sure

#    Jup = M   # Note: M here is the local M for each process
#    Jupl = Jup + 1
    msgUP = 10
    msgDN = 20
#    I2 = M + 2

    if (Me != 1):# Send bottom row to neighbor down
        msgtag = msgUP
        MPI.COMM_WORLD.Send( [U[1,:], h] , dest = NodeDN, tag = msgtag) #NodeLEFT


    if (Me != nWRs): # Receives bottom row from upper neighbor and saves as upper boundary
        msgtag = msgUP
        MPI.COMM_WORLD.Recv( [U[Mz+1,:], h] , source = NodeUP, tag = msgtag) # NodeRIGHT


    if (Me != nWRs):  # Send the top row to neighbor up
        msgtag = msgDN
        MPI.COMM_WORLD.Send( [U[Mz,:], h] , dest = NodeUP, tag = msgtag)  #NodeRight, msgtag


    if (Me != 1): # Recieves top row from lower neighbor and saves as lower boundary
        msgtag = msgDN
        MPI.COMM_WORLD.Recv( [U[0,:], h] , source = NodeDN, tag = msgtag) #NodeLEFT, msgtag

    h.Free()

    return U






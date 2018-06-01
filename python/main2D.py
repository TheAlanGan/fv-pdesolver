# -*- coding: utf-8 -*-
"""
Alan Gan - 11/29/17
M578 - 2D Parallel

main

"""
from mpi4py import MPI
import mainMR2D as mainMR
import mainWR2D as mainWR

nPROC = MPI.COMM_WORLD.Get_size()
myID = MPI.COMM_WORLD.Get_rank()

mster = 0           # rank 0 is master process
nWRs = nPROC - 1    # Number of Workers

#print "Running with ", nWRs, " workers."

if(myID == mster):
    tt0 = MPI.Wtime()            # Starts timer
    mainMR.MASTER(nWRs, mster)
    tt1 = MPI.Wtime()            # Ends timer

    print 'MR timing = ', tt1-tt0, ' sec on ', nWRs, ' workers.'
    print 'Machine: ', MPI.Get_processor_name()
    print

else:
    mainWR.WORKER(nWRs, myID)

MPI.Finalize()

# -*- coding: utf-8 -*-
"""
Alan Gan - 11/29/17
M578 - 2D Parallel

Update
"""
import numpy as np

def FLUX(RIn, ROut, ZIn, ZOut, R, Z, D, dr, dz, Mr, Mz, Fr, Fz, U, time, Me, nWRs):
### Note: U array is flipped

    """_____REAL BOUNDARY CONDITIONS_____""" # Real (not local) boundaries
    if Z[0] == ZIn: # If REAL bottom boudary, apply BC
        U[0, :]    = 0.0         # Bottom Boundary (Dirichlet)

    if Z[Mz+1] == ZOut: # If REAL top boundary, apply BC
        U[Mz+1, :] = 0.0         # Top Boundary    (Dirichlet)

    U[:, 0]    = 0.0         # Left Boundary   (Dirichlet)
    U[:, Mr+1] = np.exp(-time) * np.log(2) * np.sin(Z[:]) # Right Boundary (Dirichlet)

    """_____END OF REAL BOUNDARY CONDITIONS_____"""



    """_____REGULAR FLUXES_____""" # Updating the fluxes
#    for j in range(1, Mz+2):
#        for i in range(1, Mr+2):
#            Fr[j, i] = -D * (U[j, i] - U[j, i-1]) /(R[i] - R[i-1])  # Flux in r direction
#            Fz[j, i] = -D * (U[j, i] - U[j-1, i]) /(Z[j] - Z[j-1])  # Flux in z direction

    Fr[1:Mz+2, 1:Mr+2] = -D * (U[1:Mz+2, 1:Mr+2] - U[1:Mz+2, 0:Mr+1]) / (R[1:Mr+2] - R[0:Mr+1])
    Fz[1:Mz+2, 1:Mr+2] = -D * (U[1:Mz+2, 1:Mr+2] - U[0:Mz+1, 1:Mr+2]) / (Z[1:Mz+2, None] - Z[0:Mz+1, None])


    """_____END OF REGULAR FLUXES_____"""


    return U, Fr, Fz



def PDE(ZIn, ZOut, R, Z, Ar, Az, dr, dz, Mr, Mz, dt, Fr, Fz, U, time, Me, nWRs): #Updates U array: U(i), i=0,...,M+1

#    for j in range(1, Mz+1):
#        for i in range(1, Mr+1):

#            FluxR = Ar[j, i] * Fr[j, i] - Ar[j, i+1] * Fr[j, i+1]
#            FluxZ = Az[j, i] * Fz[j, i] - Az[j+1, i] * Fz[j+1, i]
#            U[j, i] = U[j, i] + dt * (FluxR + FluxZ) / (2*np.pi*R[i]*dz*dr)

    # Updates cell values using fluxes
    U[1:Mz+1, 1:Mr+1] = U[1:Mz+1, 1:Mr+1] + dt * ( Ar[1:Mz+1, 1:Mr+1] * Fr[1:Mz+1, 1:Mr+1] - Ar[1:Mz+1, 2:Mr+2] * Fr[1:Mz+1, 2:Mr+2] + Az[1:Mz+1, 1:Mr+1] * Fz[1:Mz+1, 1:Mr+1] - Az[2:Mz+2, 1:Mr+1] * Fz[2:Mz+2, 1:Mr+1] ) / (2*np.pi*R[None,1:Mr+1]*dz*dr)

    return U


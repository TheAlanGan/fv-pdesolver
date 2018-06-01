# -*- coding: utf-8 -*-
"""
Alan Gan - 11/29/17
M578 - 2D Parallel

Setup
"""
import numpy as np

def MESH(RIn, ROut, ZIn, ZOut, dr, dz, Mr, Mz, Me, nWRs): # Returns R and Z Mesh... Returns both Area arrays

    """_____BUILDING R MESH_____"""
    R = np.zeros(Mr + 2)  # Building the R positional Mesh
    R[0]    = RIn      # Boundary location
    R[Mr+1] = ROut     # Boundary location
    R[1]    = RIn + dr/2

    for i in range (2, Mr+1): # Standard calculation of R location
        R[i] = R[1] + (i-1) * dr
    """_____END OF R MESH_____"""



    """_____BUILDING Z MESH_____"""
    Z = np.zeros(Mz + 2)  # Building the Z positional Mesh

    Z[0]    = ZIn       # Boundary location
    Z[Mz+1] = ZOut      # Boundary location
    Z[1]    = ZIn + dz/2

    for j in range (2, Mz+1): # Standard calculation of Z location
        Z[j] = Z[1] + (j-1) * dz

    # Next part is to account for the fact that internal boundaries are differently
    # spaced from the actual boundaries (dz versus dz/2)

    if (Me == 1) and (Me != nWRs) and (Me != 0):
        Z[Mz+1] = ZOut + dz/2

    if (Me == nWRs) and (Me != 1) and (Me != 0):
        Z[0] = ZIn - dz/2

    if (Me != 1) and (Me != nWRs) and (Me != 0):
        Z[Mz+1] = ZOut + dz/2
        Z[0] = ZIn - dz/2

    """_____END OF Z MESH_____"""



    """_____BUILDING AREA MESHES_____"""
    Ar = np.zeros((Mz+2, Mr+2)) # Area of R face at R_i-1/2,j
    Az = np.zeros((Mz+2, Mr+2)) # Area of Z face at Z_i,j-1/2

    for j in range (1, Mz+2):
        for i in range (1, Mr+2):
            Ar[j, i] = 2 * np.pi * (R[i] + R[i-1])/2 * dz#(Z[j] - Z[j-1]) #dz #####(R[i] + R[i-1])/2  * dz  # Left Area at each (non-left or right bdry) node

            Az[j, i] = 2 * np.pi * R[i] * dr # Bottom Area at each (non-top or bott bdry) node)

    Ar[:, 1]    = 2 * np.pi * R[0] * dz#( ((Z[2]+Z[1])/2) - Z[0] ) #dz    # Left boundary area for left boundary nodes
    Ar[:, Mr+1] = 2 * np.pi * R[Mr+1] * dz#( Z[Mz+1] - ((Z[Mz]+Z[Mz-1])/2) ) #dz # Right boundary area for right boundary nodes

    """_____END OF AREA MESHES_____"""

    return R, Z, Ar, Az



def INIT(time, R, Z, Mr, Mz):
    """_____INITIALIZING THE ARRAYS_____"""
    U = np.zeros((Mz+2, Mr+2))   #
    Fr = np.zeros((Mz+2, Mr+2))  # The Flux in the R direction
    Fz = np.zeros((Mz+2, Mr+2))  # The Flux in the Z directions
    """_____END OF INITIALIZATION____"""


    """_____INITIAL CONDITIONS_____"""
    for j in range (0, Mz+2):
        for i in range(0, Mr+2):
            U[j, i] = np.log(R[i]) * np.sin(Z[j])   # Initial Conditions
    """_____END OF INITIAL CONDITIONS_____"""

    return U, Fr, Fz

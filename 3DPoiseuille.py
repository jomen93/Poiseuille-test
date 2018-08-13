import numpy as np
import matplotlib.pyplot as plt

# Constantes
Nx = 128; Ny = 128; Nz = 5
Nx1 = Nx+1; Ny1 = Nz+1; Nz+1
Q =19
rho0 = 1.
uxo = uyo = uzo = 0. 
g = 1e-6; dt =1.

U = uxo*np.zeros((Nx,Ny,Nz))
V = uyo*np.zeros((Nx,Ny,Nz))
W = uzo*np.zeros((Nx,Ny,Nz))

cx = np.array([0, 1,-1, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 1,-1, 1,-1, 0, 0])
cy = np.array([0, 0, 0, 1,-1, 0, 0, 1,-1, 1, 0, 1,-1,-1, 1, 0, 0, 1,-1])
cz = np.array([0, 0, 0, 0, 0, 1,-1, 0, 0, 1,-1, 1,-1, 1, 0,-1, 1,-1, 1])
w = np.array([1.0/3 ,1.0/18,1.0/18,1.0/18,1.0/18,1.0/18,
			   1.0/18,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,
			   1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36,1.0/36])


def feq(rho,U,V,W,q):

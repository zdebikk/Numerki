import numpy as np
from scipy.sparse import diags
import matplotlib.pyplot as plt
from MeshingFunctions import *
from ReferenceElement import *
from AnalyticalInfo import *

# Solution of u_t = u_{xx} in (0,1) with homogeneous Dirichlet boundary conditions

plt.rcParams['savefig.dpi'] = 300

# Problem parameters
nu = 1    # Thermal conductivity
Tf = 0.2  # Final time

## Finite element method in space
# Reference element definition
degree=1
refElement=defineRefElement1D(degree)

# Computational mesh
domain = np.array([0,1]); nOfElements=10; h = 1/nOfElements
X,T=uniform1dMesh(domain[0],domain[1],nOfElements,degree)

# System of ODE
nOfNodes = X.shape[0]; 
MprDiag = np.concatenate((2*np.ones(1),4*np.ones(nOfNodes-2),2*np.ones(1)))
MupDiag = np.ones(nOfNodes-1)
M = h/6*diags([MupDiag,MprDiag,MupDiag],[-1,0,1]).toarray()

KprDiag = np.concatenate((np.ones(1),2*np.ones(nOfNodes-2),np.ones(1)))
KupDiag = -1*np.ones(nOfNodes-1)
K = (1/h)*diags([KupDiag,KprDiag,KupDiag],[-1,0,1]).toarray()

f = np.zeros((nOfNodes,1)) 

# Homogeneous Dirichlet boundary conditions at both ends
nodesDir = np.array([0,nOfNodes-1]); 
nodesUnk = np.setdiff1d(np.arange(nOfNodes),nodesDir)
valDir = np.zeros((len(nodesDir),1))
dtvalDir = np.zeros((len(nodesDir),1))

# System reduction
f_red = f[nodesUnk] - K[nodesUnk,:][:,nodesDir]@valDir - M[nodesUnk,:][:,nodesDir]@dtvalDir
K_red = K[nodesUnk,:][:,nodesUnk]
M_red = M[nodesUnk,:][:,nodesUnk]

## Euler explicit time intergration
# Time steps
timeStep = 0.001
#c
timeStep = h**2/(6*nu)*0.5 #last 0.5 is to ensure3 its smaller

timeIterationsExp = round(Tf/timeStep+1); timeStep = Tf/timeIterationsExp
timeVec = np.linspace(0,Tf,timeIterationsExp); 
print('Time steps',timeIterationsExp)
print('Euler explicit: nu*At/h^2=', nu*timeStep*nOfElements**2)
# Allocation of the solution at each time-step
u_SolExp = np.zeros((nOfNodes,timeIterationsExp))
# Initial condition
inCond = initial_Cond(X); inCond.shape = (nOfNodes)
u_SolExp[:,0] = inCond

# Euler explicit time integration (from exercise)
for k in np.arange(timeIterationsExp-1):
    tk = timeVec[k]
    tk1 = timeVec[k+1]
    aux1 = -nu*K_red@u_SolExp[1:nOfNodes-1,k]; aux1.shape = (nOfNodes-2,1)
    f_2 = aux1 + f_red
    aux = np.linalg.solve(M_red,f_2); aux.shape = nOfNodes-2
    u_SolExp[1:nOfNodes-1,k+1] = u_SolExp[1:nOfNodes-1,k] + timeStep*aux

#Backward Euler
u_SolExp_BE = np.zeros((nOfNodes,timeIterationsExp))
u_SolExp_BE[:,0] = inCond

# Euler backward time integration
for k in np.arange(timeIterationsExp-1):
    tk = timeVec[k]
    tk1 = timeVec[k+1]
    aux1 = M_red@u_SolExp_BE[1:nOfNodes-1,k]; aux1.shape = (nOfNodes-2,1)
    f_2 = aux1 + f_red*timeStep
    BEMATRIX = M_red + timeStep*nu*K_red
    aux = np.linalg.solve(BEMATRIX,f_2); aux.shape = nOfNodes-2
    u_SolExp_BE[1:nOfNodes-1,k+1] = aux

#a
print(u_SolExp[:,0])
print(u_SolExp[5,0])
print(u_SolExp[5,20])
print(u_SolExp[5,timeIterationsExp-1]) #solution at final time

# num = 15
# for k in np.arange(round(timeIterationsExp/num)):
#     plt.plot(X,u_SolExp[:,num*k],'-o')
# plt.show()

plt.imshow(u_SolExp ,aspect='auto', cmap='viridis')
plt.colorbar()
plt.title('Euler Explicit', fontsize=14, pad=12)
plt.show()

plt.imshow(u_SolExp, aspect='auto', interpolation='bicubic', cmap='viridis')
plt.colorbar()
plt.title('Euler Explicit (Bicubic Interpolation)', fontsize=14, pad=12)
plt.show()

plt.imshow(u_SolExp_BE ,aspect='auto', cmap='viridis')
plt.colorbar()
plt.title('Backward Euler', fontsize=14, pad=12)
plt.show()

plt.imshow(u_SolExp - u_SolExp_BE ,aspect='auto', cmap='viridis')
plt.colorbar()
plt.title('Difference Euler Explicit - Backward Euler', fontsize=14, pad=12)
plt.show()




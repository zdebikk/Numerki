import numpy as np
from scipy.sparse import diags
import matplotlib.pyplot as plt
from MeshingFunctions import *
from ReferenceElement import *
from AnalyticalInfo import *
import sys
from pathlib import Path
# Add 'mainfolder' (2 levels up from this file) to Python's search path
root_dir = Path(__file__).resolve().parent.parent
if str(root_dir) not in sys.path:
    sys.path.append(str(root_dir))
# Now standard import works
from Sergi.MeshingFunctions import UniformRectangleMesh
from Sergi.ReferenceElement import defineRefElement
# from Sergi.main_Laplace2D_forAlba import laplace2D_solver
from Sergi.PlottingFunctions import plotMesh
from Sergi.SystemComputation import computeSystemLaplace



# Solution of u_t = u_{xx} in (0,1) with homogeneous Dirichlet boundary conditions

plt.rcParams['savefig.dpi'] = 300

# Problem parameters
nu = 1    # Thermal conductivity
Tf = 0.2  # Final time

## Finite element method in space
# Reference element definition\
degree = 1
refElement=defineRefElement(elementType = 0, degree = degree)

RESOL_OF_SPACE_DOMAIN = 10 # Number of elements in the space domain
# Computational mesh
domain = np.array([0,5,0,5])
nx = RESOL_OF_SPACE_DOMAIN; ny = nx; h = max(domain[1]-domain[0],domain[3]-domain[2])/nx;
X,T = UniformRectangleMesh(domain,nx,ny,refElement)

plotMesh(X,T,refElement); plt.show()

print(X)

print(T)

[K,f] = computeSystemLaplace(X,T,refElement)


# Dirichlet boundary conditions
x1 = domain[0]; x2 = domain[1]
y1 = domain[2]; y2 = domain[3]
tol = 1e-8
nodes_Left = np.where(abs(X[:,0]-x1)<tol)[0]
nodes_Right = np.where(abs(X[:,0]-x2)<tol)[0]
nodes_Bottom = np.where(abs(X[:,1]-y1)<tol)[0]
nodes_Top = np.where(abs(X[:,1]-y2)<tol)[0]
nodesDir = np.unique(np.block([nodes_Left,nodes_Right,nodes_Bottom,nodes_Top]))
valDir = boundary_Cg(X[nodesDir,:])[0]
valDir.shape = (len(nodesDir),1)

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




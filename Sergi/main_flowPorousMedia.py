import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from MeshingFunctions import *
from PlottingFunctions import *
from ReferenceElement import * 
from SystemComputation import *

plt.rcParams['savefig.dpi'] = 300


do_plot = 1

# Reference element
elementType = 0 # 0: quadrilaterals, 1 = triangles
degree = 1
refElement = defineRefElement(elementType,degree)

# Computational mesh
fineness = 1 #from 0 to 3
wallDistance = 10 # 5 or 10
X,T = ExcavationMeshes(degree,elementType,fineness,wallDistance)
if do_plot == 1:
    plotMesh(X,T,refElement); plt.show()

# FE system assembly
[K,f] = computeSystemPorousMedia(X,T,refElement)
#[K,f] = computeSystemLaplaceSparse(X,T,referenceElement)
plt.spy


# Dirichlet boundary conditions
tol = 1e-8
nodes_ex1 = np.where((abs(X[:,1]-1.5)<tol) & (abs(X[:,0]-0.75)<tol+0.75))[0]
nodes_ex2 = np.where((abs(X[:,1]-3)<tol) & (X[:,0]>1.75-tol))[0]
nodesDir = np.unique(np.block([nodes_ex1,nodes_ex2]))
valDir = exactSol2(X[nodesDir,:])[0]
valDir.shape = (len(nodesDir),1)

# System reduction and solution 
u = findSolution_SystemReduction(K,f,nodesDir,valDir)

contourfPlot(u,X,T) 
plt.gcf().canvas.manager.set_window_title('Contour Plot - Numerical Solution u')
plt.savefig('Figure_ns.png', dpi=300, bbox_inches='tight')
plt.show()

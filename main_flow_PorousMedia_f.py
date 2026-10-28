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

mesh_sizes = []
flow_values = []
flow_values2 = []
wallDistance = 10 # 10m is chosen based on observation 2(d)

for fineness in range(4):
    X, T = ExcavationMeshes(degree, elementType, fineness, wallDistance)
    
    nElements = T.shape[0]
    h_approx = 1.0 / np.sqrt(nElements) 
    
    [K, f] = computeSystemPorousMedia(X, T, refElement)
    
    tol = 1e-8
    nodes_ex1 = np.where((abs(X[:,1]-1.5)<tol) & (abs(X[:,0]-0.75)<tol+0.75))[0]
    nodes_ex2 = np.where((abs(X[:,1]-3)<tol) & (X[:,0]>1.75-tol))[0]
    nodesDir = np.unique(np.block([nodes_ex1, nodes_ex2]))
    valDir = exactSol2(X[nodesDir,:])[0].reshape(-1, 1)
    
    u_lag, lambda_lag = findSolution_LagrangeMultipliers(K, f, nodesDir, valDir)
    
    total_flow = np.abs(np.sum(lambda_lag))
    
    mesh_sizes.append(h_approx)
    flow_values.append(total_flow)
    print(fineness)


plt.figure(figsize=(8, 5))
plt.plot(mesh_sizes, flow_values, 'o-')
plt.title('Influence of Mesh Size on Excavation Surface Flow')
plt.xlabel('Approximate Mesh Size $h \propto 1/\sqrt{N_{elements}}$')
plt.ylabel('Total Flow $Q$ ($m^3/s$)')
plt.grid(True, linestyle='--')
plt.gca().invert_xaxis()
plt.show()

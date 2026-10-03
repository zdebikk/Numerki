import numpy as np
from scipy.sparse import diags
from analyticalInfo import * 
from  PlottingFunctions import *
from ReferenceElement import *
from ErrorComputation import *
from MeshingFunctions import *
from SystemComputation import *
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
# This library enhances some plot options
from mpl_toolkits.mplot3d import Axes3D
# Plots are shown on a new window, where we can move and zoom
#matplotlib auto  

do_plot = 1

# Problem definition 
domain = np.array([0,1,0,1])

def boundary_func(x, y):
    return np.sin(15 * x+1) * np.cos(9 * y) + 2*x*y + 4*2**((-((x-0.5)**2 + (y-0.5)**2))/1e-2)

elementType = 0 # 1 = Triangles, 0 = Quadrilaterals
degree = 1
referenceElement = defineRefElement(elementType, degree)

# Creation of the mesh
nx = 32; ny = nx; h = max(domain[1]-domain[0],domain[3]-domain[2])/nx; # Number of elements in each direction and element size 
print('Number of elements',np.array([nx,ny]))
X,T = UniformRectangleMesh(domain,nx,ny,referenceElement)
# if do_plot == 1:
    #plotMesh(X,T,referenceElement); plt.show()

# FE system assembly
[K,f] = computeSystemLaplace(X,T,referenceElement)
#[K,f] = computeSystemLaplaceSparse(X,T,referenceElement)
#plt.spy(K); plt.show()

f = np.zeros_like(f)

# Dirichlet boundary conditions
x1 = domain[0]; x2 = domain[1]
y1 = domain[2]; y2 = domain[3]
tol = 1e-8
tol2 = 1e-2
nodes_Left = np.where(abs(X[:,0]-x1)<tol)[0]
nodes_Right = np.where(abs(X[:,0]-x2)<tol)[0]
nodes_Bottom = np.where(abs(X[:,1]-y1)<tol)[0]
nodes_Top = np.where(abs(X[:,1]-y2)<tol)[0]
nodes_Middle = np.where((abs(X[:,0]-(x1+x2)/2) < tol2) & (abs(X[:,1]-(y1+y2)/2) < tol2))[0]
nodesDir = np.unique(np.block([nodes_Left,nodes_Right,nodes_Bottom,nodes_Top,nodes_Middle]))
x_coords = X[nodesDir, 0]
y_coords = X[nodesDir, 1]
valDir = boundary_func(x_coords, y_coords).reshape((len(nodesDir), 1))

# System reduction and solution 
u = findSolution_SystemReduction(K,f,nodesDir,valDir)
#u = findSolutionSparse_SystemReduction(K,f,nodesDir,valDir)





nOfNodes = np.shape(X)[0]
if do_plot==1:
    # Contour plot
    contourfPlot(u,X,T) 
    plt.gcf().canvas.manager.set_window_title('Contour Plot - Numerical Solution u')
    plt.show()
    # Surface plot
    surfPlot(u,X,T)
    plt.gcf().canvas.manager.set_window_title('Surface Plot - Numerical Solution u')
    plt.show()




import numpy as np
import scipy.io as sio
import matplotlib.pyplot as plt
from MeshingFunctions import *
from PlottingFunctions import *
from ReferenceElement import * 
from SystemComputation import *

do_plot = 1

# Reference element
elementType = 0 # 0: quadrilaterals, 1 = triangles
degree = 2
refElement = defineRefElement(elementType,degree)

# Computational mesh
fineness = 2 #from 0 to 3
wallDistance = 5 # 5 or 10
X,T = ExcavationMeshes(degree,elementType,fineness,wallDistance)
if do_plot == 1:
    plotMesh(X,T,refElement); plt.show()
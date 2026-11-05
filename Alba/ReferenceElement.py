import numpy as np
import sys

def defineRefElement1D(degree):
    class refElement:
        pass
    refElement.degree = degree
    if degree == 1:
        refElement.nOfElementNodes = 2
        refElement.nodesCoord = np.array([-1,1])
        a = np.sqrt(1/3)
        xi = np.array([[-a],[a]])
        refElement.integrationPoints = xi
        refElement.integrationWeights = np.array([1,1])
        refElement.nIP = 2
        refElement.N = np.block([(1-xi)/2,(1+xi)/2])
        refElement.Nxi = np.array([[-1,1],[-1,1]])
    else:
        print('Not implemented element')
        sys.exit()   
        
    return refElement
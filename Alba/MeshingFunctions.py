import numpy as np
import sys

def uniform1dMesh(a,b,nOfElements,degree): 
    X = np.linspace(a,b,nOfElements+1)
    X.shape = (nOfElements+1,1)
    if degree == 1:
        aux=np.arange(0,nOfElements)
        aux.shape=(nOfElements,1)
        T=np.block([aux,aux+1])
    else:
        print('Not implemented element')
        sys.exit()   
                 
    return X,T    
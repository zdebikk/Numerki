import numpy as np
from AnalyticalInfo import *

def computeSystem_HeatEqn1D(X,T):
    nOfNodes = X.shape[0]
    nOfElements = T.shape[0]; nOfElementNodes = T.shape[1]
    K=np.zeros((nOfNodes,nOfNodes)) #full matrix, sparse storage in 2D...
    M=np.zeros((nOfNodes,nOfNodes))
    f=np.zeros((nOfNodes,1))
    #Loop in elements for matrix computation
    for e in np.arange(nOfElements):
        he = X[e+1]-X[e]
        #Elemental matrix for 1D linear element
        Ke = np.array([[1,-1],[-1,1]])*(1/he)
        Me = np.array([[2,1],[1,2]])*(he/6)
        #Matrix assembly
        for i in np.arange(nOfElementNodes):
            for j in np.arange(nOfElementNodes):
                K[T[e,i],T[e,j]]=K[T[e,i],T[e,j]]+Ke[i,j]
                M[T[e,i],T[e,j]]=M[T[e,i],T[e,j]]+Me[i,j]  
        MLumped = np.zeros((nOfNodes,nOfNodes))
        for i in np.arange(nOfNodes):
            MLumped[i,i] = sum(M[i,:]) 
    return K,M,MLumped,f
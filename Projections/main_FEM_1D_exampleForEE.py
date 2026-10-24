import numpy as np
from scipy.sparse import diags
import matplotlib.pyplot as plt
from femFunctions1D import *

#Definition of source term s=-u_xx and Dirichlet values at both end points
def sourceTerm(x):
    return 1/5+20*(-200*x**2 +280*x -97)*np.exp(-(10*x-7)**2)
ua=0
ub=0

degree=1
referenceElement=defineRefElement1D(degree)

#uniform mesh
nOfElements=10
X,T=uniform1dMesh(0,1,nOfElements,degree)

nOfNodes=X.shape[0]
nOfElementNodes=T.shape[1]
K=np.zeros((nOfNodes,nOfNodes)) #full matrix, sparse storage in 2D...
f=np.zeros((nOfNodes,1))
#Loop in elements for matrix computation
for e in np.arange(nOfElements):
    he=X[e+1]-X[e]
    #Elemental matrix for 1D linear element
    Ke=np.array([[1,-1],[-1,1]])*(1/he)
    #Matrix assembly
    for i in np.arange(nOfElementNodes):
        for j in np.arange(nOfElementNodes):
            K[T[e,i],T[e,j]]=K[T[e,i],T[e,j]]+Ke[i,j]
#Loop in elements for rhs vector computation
xi=referenceElement.integrationPoints
w=referenceElement.integrationWeights
nOfIntegrationPoints=xi.shape[0]
for e in np.arange(nOfElements):
    he=X[e+1]-X[e]
    #Elemental vector with numerical integration
    fe=np.zeros((1,nOfElementNodes))
    for k in np.arange(nOfIntegrationPoints):
        Nk=referenceElement.N[k,:]
        xk=(X[e]+X[e+1])/2+he/2*xi[k] #linear transformation
        fe=fe+Nk*(sourceTerm(xk)*(he/2)*w[k])
    fe.shape=(nOfElementNodes,1)
    #Vector assembly
    for i in np.arange(nOfElementNodes):
        f[T[e,i]]= f[T[e,i]] + fe[i]
    
#Dirichlet boundary conditions at both ends -> system reduction
aux=K[1:-1,0]*ua+K[1:-1,-1]*ub
aux.shape=(nOfNodes-2,1)
freduced=f[1:-1]-aux
Kreduced=K[1:-1,1:-1]

#System solution
solution=np.linalg.solve(Kreduced,freduced)

u=np.zeros((nOfNodes,1))
u[0]=ua
u[1:-1]=solution
u[-1]=ub

# print(np.linalg.solve(np.array([[1, 1], [0, 1]]), np.array([2, 3])))

#

#Postprocess
plt.plot(X,u,'-o')
plt.xlabel('x')
plt.ylabel('u')
plt.show()

#Derivative in each element
dudx=(u[1:]-u[:-1])/(X[1:]-X[:-1])
plt.xlabel('x')
plt.ylabel('du/dx')
plotElementalConstants(X,T,dudx)

#Computation of the L2 gradient smoothing (L2 projection of du/dx in the FE space)

#Example of refinement
elementsToRefine=np.array([1,4])
Xnew,Tnew=refineMesh1DLinear(X,T,elementsToRefine)
print('Xnew=',Xnew)






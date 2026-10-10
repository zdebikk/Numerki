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

remote = True

do_plot = 0

# Problem definition 
domain = np.array([0,1,0,1])

hsarray = list(range(2,15))

if remote:
    hsarray = []
    for exp in range(16,18):
        hsarray.append(int(1.3**(exp)))


hsdata = [[_] for _ in hsarray]

print("hsarray:", hsarray)




# print("KWADRATY 1 STOPNIA")
# increment = -1
# for hs in hsarray:
#     increment += 1
#     # Reference element and computational mesh
#     elementType = 0 # 1 = Triangles, 0 = Quadrilaterals
#     degree = 1
#     referenceElement = defineRefElement(elementType, degree)

#     # Creation of the mesh
#     nx = hs; ny = nx; h = max(domain[1]-domain[0],domain[3]-domain[2])/nx; # Number of elements in each direction and element size 
#     print('Number of elements',np.array([nx,ny]))
#     X,T = UniformRectangleMesh(domain,nx,ny,referenceElement)
#     if do_plot == 1:
#         plotMesh(X,T,referenceElement); plt.show()

#     # FE system assembly
#     [K,f] = computeSystemLaplace(X,T,referenceElement)
#     #[K,f] = computeSystemLaplaceSparse(X,T,referenceElement)
#     #plt.spy(K); plt.show()

#     # Dirichlet boundary conditions
#     x1 = domain[0]; x2 = domain[1]
#     y1 = domain[2]; y2 = domain[3]
#     tol = 1e-8
#     nodes_Left = np.where(abs(X[:,0]-x1)<tol)[0]
#     nodes_Right = np.where(abs(X[:,0]-x2)<tol)[0]
#     nodes_Bottom = np.where(abs(X[:,1]-y1)<tol)[0]
#     nodes_Top = np.where(abs(X[:,1]-y2)<tol)[0]
#     nodesDir = np.unique(np.block([nodes_Left,nodes_Right,nodes_Bottom,nodes_Top]))
#     valDir = exactSol(X[nodesDir,:])[0];valDir.shape = (len(nodesDir),1)

#     # System reduction and solution 
#     u = findSolution_SystemReduction(K,f,nodesDir,valDir)
#     #u = findSolutionSparse_SystemReduction(K,f,nodesDir,valDir)

#     # L2 error computation
#     L2Error = computeL2Error(u,X,T,referenceElement)
#     print('L2 error: ', L2Error)

#     # H1 error computation
#     H1Error = computeH1Error(u,X,T,referenceElement)
#     print('H1 error: ', H1Error)

#     hsdata[increment].append(L2Error[0])
#     hsdata[increment].append(H1Error[0])


# print("TROJKATY 1 STOPNIA")
# increment = -1
# for hs in hsarray:
#     increment += 1
#     # Reference element and computational mesh
#     elementType = 1 # 1 = Triangles, 0 = Quadrilaterals
#     degree = 1
#     referenceElement = defineRefElement(elementType, degree)

#     # Creation of the mesh
#     nx = hs; ny = nx; h = max(domain[1]-domain[0],domain[3]-domain[2])/nx; # Number of elements in each direction and element size 
#     print('Number of elements',np.array([nx,ny]))
#     X,T = UniformRectangleMesh(domain,nx,ny,referenceElement)
#     if do_plot == 1:
#         plotMesh(X,T,referenceElement); plt.show()

#     # FE system assembly
#     [K,f] = computeSystemLaplace(X,T,referenceElement)
#     #[K,f] = computeSystemLaplaceSparse(X,T,referenceElement)
#     #plt.spy(K); plt.show()

#     # Dirichlet boundary conditions
#     x1 = domain[0]; x2 = domain[1]
#     y1 = domain[2]; y2 = domain[3]
#     tol = 1e-8
#     nodes_Left = np.where(abs(X[:,0]-x1)<tol)[0]
#     nodes_Right = np.where(abs(X[:,0]-x2)<tol)[0]
#     nodes_Bottom = np.where(abs(X[:,1]-y1)<tol)[0]
#     nodes_Top = np.where(abs(X[:,1]-y2)<tol)[0]
#     nodesDir = np.unique(np.block([nodes_Left,nodes_Right,nodes_Bottom,nodes_Top]))
#     valDir = exactSol(X[nodesDir,:])[0];valDir.shape = (len(nodesDir),1)

#     # System reduction and solution 
#     u = findSolution_SystemReduction(K,f,nodesDir,valDir)
#     #u = findSolutionSparse_SystemReduction(K,f,nodesDir,valDir)

#     # L2 error computation
#     L2Error = computeL2Error(u,X,T,referenceElement)
#     print('L2 error: ', L2Error)

#     # H1 error computation
#     H1Error = computeH1Error(u,X,T,referenceElement)
#     print('H1 error: ', H1Error)

#     hsdata[increment].append(L2Error[0])
#     hsdata[increment].append(H1Error[0])




print("TROJKATY 2 STOPNIA")
increment = -1
for hs in hsarray:
    increment += 1
    # Reference element and computational mesh
    elementType = 1 # 1 = Triangles, 0 = Quadrilaterals
    degree = 2
    referenceElement = defineRefElement(elementType, degree)

    # Creation of the mesh
    nx = hs; ny = nx; h = max(domain[1]-domain[0],domain[3]-domain[2])/nx; # Number of elements in each direction and element size 
    print('Number of elements',np.array([nx,ny]))
    X,T = UniformRectangleMesh(domain,nx,ny,referenceElement)
    if do_plot == 1:
        plotMesh(X,T,referenceElement); plt.show()

    # FE system assembly
    [K,f] = computeSystemLaplace(X,T,referenceElement)
    #[K,f] = computeSystemLaplaceSparse(X,T,referenceElement)
    #plt.spy(K); plt.show()

    # Dirichlet boundary conditions
    x1 = domain[0]; x2 = domain[1]
    y1 = domain[2]; y2 = domain[3]
    tol = 1e-8
    nodes_Left = np.where(abs(X[:,0]-x1)<tol)[0]
    nodes_Right = np.where(abs(X[:,0]-x2)<tol)[0]
    nodes_Bottom = np.where(abs(X[:,1]-y1)<tol)[0]
    nodes_Top = np.where(abs(X[:,1]-y2)<tol)[0]
    nodesDir = np.unique(np.block([nodes_Left,nodes_Right,nodes_Bottom,nodes_Top]))
    valDir = exactSol(X[nodesDir,:])[0];valDir.shape = (len(nodesDir),1)

    # System reduction and solution 
    u = findSolution_SystemReduction(K,f,nodesDir,valDir)
    #u = findSolutionSparse_SystemReduction(K,f,nodesDir,valDir)

    # L2 error computation
    L2Error = computeL2Error(u,X,T,referenceElement)
    print('L2 error: ', L2Error)

    # H1 error computation
    H1Error = computeH1Error(u,X,T,referenceElement)
    print('H1 error: ', H1Error)

    hsdata[increment].append(L2Error[0])
    hsdata[increment].append(H1Error[0])




# print("KWADRATY 2 STOPNIA")
# increment = -1
# for hs in hsarray:
#     increment += 1
#     # Reference element and computational mesh
#     elementType = 0 # 1 = Triangles, 0 = Quadrilaterals
#     degree = 2
#     referenceElement = defineRefElement(elementType, degree)

#     # Creation of the mesh
#     nx = hs; ny = nx; h = max(domain[1]-domain[0],domain[3]-domain[2])/nx; # Number of elements in each direction and element size 
#     print('Number of elements',np.array([nx,ny]))
#     X,T = UniformRectangleMesh(domain,nx,ny,referenceElement)
#     if do_plot == 1:
#         plotMesh(X,T,referenceElement); plt.show()

#     # FE system assembly
#     [K,f] = computeSystemLaplace(X,T,referenceElement)
#     #[K,f] = computeSystemLaplaceSparse(X,T,referenceElement)
#     #plt.spy(K); plt.show()

#     # Dirichlet boundary conditions
#     x1 = domain[0]; x2 = domain[1]
#     y1 = domain[2]; y2 = domain[3]
#     tol = 1e-8
#     nodes_Left = np.where(abs(X[:,0]-x1)<tol)[0]
#     nodes_Right = np.where(abs(X[:,0]-x2)<tol)[0]
#     nodes_Bottom = np.where(abs(X[:,1]-y1)<tol)[0]
#     nodes_Top = np.where(abs(X[:,1]-y2)<tol)[0]
#     nodesDir = np.unique(np.block([nodes_Left,nodes_Right,nodes_Bottom,nodes_Top]))
#     valDir = exactSol(X[nodesDir,:])[0];valDir.shape = (len(nodesDir),1)

#     # System reduction and solution 
#     u = findSolution_SystemReduction(K,f,nodesDir,valDir)
#     #u = findSolutionSparse_SystemReduction(K,f,nodesDir,valDir)

#     # L2 error computation
#     L2Error = computeL2Error(u,X,T,referenceElement)
#     print('L2 error: ', L2Error)

#     # H1 error computation
#     H1Error = computeH1Error(u,X,T,referenceElement)
#     print('H1 error: ', H1Error)

#     hsdata[increment].append(L2Error[0])
#     hsdata[increment].append(H1Error[0])







# nOfNodes = np.shape(X)[0]
# if do_plot==1:
#     # Contour plot
#     contourfPlot(u,X,T) 
#     plt.show()
#     u_analytic = exactSol(X)[0].reshape(nOfNodes,1)
#     contourfPlot(u_analytic,X,T)
#     plt.show()
#     # Surface plot
#     surfPlot(u,X,T)
#     plt.show()
#     surfPlot(u_analytic,X,T)
#     plt.show()

print('hsdata:', hsdata)


if not remote:

    # hsdata = [[2, np.float64(0.6410725737160661), np.float64(3.93918016803861), np.float64(0.6323877383818364), np.float64(5.001707267371442)], [3, np.float64(0.7199909266191815), np.float64(3.4780028565343324), np.float64(0.2578187839505659), np.float64(5.310846236668786)], [4, np.float64(0.5024967946510979), np.float64(4.320852553232804), np.float64(0.2542965394764775), np.float64(4.394112872316398)], [5, np.float64(0.3095720577436989), np.float64(4.4367912894434385), np.float64(0.17571760055679045), np.float64(3.174680219840116)], [6, np.float64(0.19297664863230088), np.float64(4.283913466052776), np.float64(0.11882294088051241), np.float64(2.4068673882993443)], [7, np.float64(0.13565194882808013), np.float64(4.079977010394338), np.float64(0.08336093985347584), np.float64(1.9071329749194337)], [8, np.float64(0.1099642815532222), np.float64(3.8670048079598796), np.float64(0.060612068659598194), np.float64(1.5497622397359576)], [9, np.float64(0.09689887147741132), np.float64(3.656864010877931), np.float64(0.04530249528668447), np.float64(1.2816760422706583)], [10, np.float64(0.08796210778106875), np.float64(3.4554342903602784), np.float64(0.034619520072264445), np.float64(1.0755627645840486)], [11, np.float64(0.08038337613188866), np.float64(3.265747404178508), np.float64(0.026963675308169475), np.float64(0.9142084842244419)], [12, np.float64(0.07345542170735586), np.float64(3.089049032306106), np.float64(0.021357102761101678), np.float64(0.785856568713227)], [13, np.float64(0.06705632554167976),np.float64(2.9255104311484086), np.float64(0.017172694486810816), np.float64(0.6822555796060302)], [14, np.float64(0.06118782183673903), np.float64(2.774696520693026), np.float64(0.013995725733391853), np.float64(0.5975207338260018)], [13, np.float64(0.06705632554167976), np.float64(2.9255104311484086), np.float64(0.017172694486810816), np.float64(0.6822555796060302)], [17, np.float64(0.046714310059389204), np.float64(2.3904888895473855), np.float64(0.008109753279989957), np.float64(0.41918687806177113)], [23, np.float64(0.02871115775505216), np.float64(1.85109075188008), np.float64(0.003400316709395511), np.float64(0.23740837565800849)], [30, np.float64(0.01791909359722215), np.float64(1.4543708069969312), np.float64(0.0015640126027480022), np.float64(0.1422874123935962)], [39, np.float64(0.010985167244123312), np.float64(1.1353454143717616), np.float64(0.0007211333009797918), np.float64(0.08520130710859407)]]

    hsdata =   [[2, np.float64(0.64107257), np.float64(3.93918017), np.float64(0.63238774), np.float64(5.00170727)], [3, np.float64(0.71999093), np.float64(3.47800286), np.float64(0.25781878), np.float64(5.31084624)], [4, np.float64(0.50249679), np.float64(4.32085255), np.float64(0.25429654), np.float64(4.39411287)], [5, np.float64(0.30957206), np.float64(4.43679129), np.float64(0.1757176), np.float64(3.17468022)], [6, np.float64(0.19297665), np.float64(4.28391347), np.float64(0.11882294), np.float64(2.40686739)], [7, np.float64(0.13565195), np.float64(4.07997701), np.float64(0.08336094), np.float64(1.90713297)], [8, np.float64(0.10996428), np.float64(3.86700481), np.float64(0.06061207), np.float64(1.54976224)], [9, np.float64(0.09689887), np.float64(3.65686401), np.float64(0.0453025), np.float64(1.28167604)], [10, np.float64(0.08796211), np.float64(3.45543429), np.float64(0.03461952), np.float64(1.07556276)], [11, np.float64(0.08038338), np.float64(3.2657474), np.float64(0.02696368), np.float64(0.91420848)], [12, np.float64(0.07345542), np.float64(3.08904903), np.float64(0.0213571), np.float64(0.78585657)], [13, np.float64(0.06705633), np.float64(2.92551043), np.float64(0.01717269), np.float64(0.68225558)], [14, np.float64(0.06118782), np.float64(2.77469652), np.float64(0.01399573), np.float64(0.59752073)], [13, np.float64(0.06705633), np.float64(2.92551043), np.float64(0.01717269), np.float64(0.68225558)], [17, np.float64(0.04671431), np.float64(2.39048889), np.float64(0.00810975), np.float64(0.41918688)], [23, np.float64(0.02871116), np.float64(1.85109075), np.float64(0.00340032), np.float64(0.23740838)], [30, np.float64(0.01791909), np.float64(1.45437081), np.float64(0.00156401), np.float64(0.14228741)], [39, np.float64(0.01098517), np.float64(1.13534541), np.float64(0.00072113), np.float64(0.08520131)], [51, np.float64(0.00656346), np.float64(0.87608725), np.float64(0.00032516), np.float64(0.05019297)]]


    x, y3, y4, y5, y6 = zip(*hsdata)

    def fun(x):
        fun.name = 'log'
        return np.log(x)
        fun.name = 'id'
        return x

    def fun2(x):
        fun2.name = 'log'
        return np.log(x)
        fun2.name = 'id'
        return x

    plt.figure(figsize=(16, 10))
    # plt.plot(fun2(x), fun(y1), label='L2 S1', marker='s', linestyle='-_',color='crimson')
    # plt.plot(fun2(x), fun(y2), label='H1 S1', marker='s', linestyle='-',color='crimson')
    plt.plot(fun2(x), fun(y3), label='L2 T1', marker='^', linestyle='--',color='royalblue')
    plt.plot(fun2(x), fun(y4), label='H1 T1', marker='^', linestyle='-',color='royalblue')
    plt.plot(fun2(x), fun(y5), label='L2 T2', marker='v', linestyle='--',color='orange')
    plt.plot(fun2(x), fun(y6), label='H1 T2', marker='v', linestyle='-',color='orange')
    # plt.plot(fun2(x), fun(y7), label='L2 S2', marker='D', linestyle='-_',color='forestgreen')
    # plt.plot(fun2(x), fun(y8), label='H1 S2', marker='D', linestyle='-',color='forestgreen')


    plt.xlabel(f'Mesh Size {fun2.name}(1/h)')
    plt.ylabel(f'{fun.name}(Errors)')
    plt.title('Errors vs. Mesh Size')
    plt.legend()
    plt.grid(True)

    # Wyświetlenie
    plt.show()
import numpy as np
def initial_Cond(x):
     return 1 - 2*np.abs(x-0.5)

def initial_CondConv(x):
     pi = np.pi
     return np.sin(pi*x)

def exactSol(x,t):
     pi = np.pi
     return np.sin(pi*x)*np.exp((-pi**2)*t)
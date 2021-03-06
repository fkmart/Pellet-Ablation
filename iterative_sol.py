import numpy as np
from gen_var import dr
from numba import jit,float64, njit
#import time 

@njit
def SOR(A,x,f,r):
    rel_tol = 1e-6
    iteration = 0 
    l = len(x)
    h = r[1] - r[0]
    omega = 2.0 / (1.0 + np.sin(np.pi * h))
    res = np.ones(len(x))
    while (np.any(res[1:-1] > np.abs(np.multiply(rel_tol,x[1:-1])))):
        for i in range(1, l-1):
            s = np.dot(A[i,:], x[:])
            xnew = f[i]*h**2  - (s - A[i,i]*x[i])
            xold = x[i]
            x[i] = (1.0 - omega)*x[i] + xnew*omega/A[i,i]
            res[i] = np.abs(xold - x[i])
        iteration += 1 
    print(iteration)
    return x 

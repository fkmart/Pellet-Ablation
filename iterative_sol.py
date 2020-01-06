import numpy as np
from gen_var import dr
from numba import jit,float64
#import time 

#@jit(nopython = True)
def SOR(A,x,f,r):
    rel_tol = 1e-8
    iteration = 0 
    omega = 1.00
    l = len(x)
    h = r[-1]/l
    h = r[1] - r[0]
    omega = 2.0 / (1.0 + np.sin(np.pi * h))

    res = np.ones(len(x))
    while (np.any(res[1:-1] > np.abs(np.multiply(rel_tol,x[1:-1])))):
        for i in range(1, l-1):
            s = np.dot(A[i,:], x[:])
            xnew = f[i]*h**2  - (s - A[i,i]*x[i]) 
            xold = x[i]
            x[i] = (1 - omega)*x[i] + xnew*omega/A[i,i]
            res[i] = np.abs(xold - x[i])
        iteration += 1
        if (iteration%10000 == 0):
            #print("\r"+ str(iteration),end='')
            print(iteration)
        else:
            pass 
    #print("\r"+ str(iteration))
    print(iteration)
    """Determining answer via matrix inversion for comparison
    -can choose to remove this for larger matrices when comparison isn't needed."""
    """
    inv = np.linalg.inv(A)
    b = np.zeros(l)
    b[0] = x[0]
    b[-1] = x[-1]
    ans = np.dot(inv,-f*h*h - b)"""
    return x 

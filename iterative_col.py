import numpy as np
from gen_var import dr
from numba import jit,float64
import time 

@jit(nopython = True)
def SOR(A,x,f,r):
    rel_tol = 1e-1
    iteration = 0 
    omega = 1.00
    r0 = 1e-3
    h = 1.0
    #omega = 2.0 / (1.0 + np.sin(np.pi * h))
    l = len(x)
    h = 1.0/l
    h = r[1] - r[0]

    res = np.ones(len(x))
    b = np.zeros(len(x))
    while (np.any(res > np.abs(np.multiply(rel_tol,x)))):
        for i in range(1, l-1):
            s = np.dot(A[i,:], x[:])
            xnew = -f[i]*h**2  - (s - A[i,i]*x[i]) - 2.0
            xold = x[i]
            x[i] = xnew*omega/A[i,i]
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
    
    inv = np.linalg.inv(A)
    ans = np.dot(inv,-f*h*h - b)
    return x, ans 

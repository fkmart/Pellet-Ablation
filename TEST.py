import numpy as np
from gen_var import dr
from numba import jit,float64
import time 

A = np.asarray([[5,0,-2], [3,5,1], [0, -3 ,4]])
f = [7,2,-4]
y = np.ones(3)


"""def th(A,y,f):
 rel_tol = 1e-6 
 iteration = 0
 omega = 1.0
 res = np.ones(len(y))
 
 while(np.any(res > np.abs(np.multiply(rel_tol,y)))):
        for i in range(0,len(y)):
            s = np.dot(A[i,:],y[:])
            # sum term for in SOR residual
            Ynew = f[i] - (s - A[i,i]*y[i]) #residual, function - boundary conditions + summed terms from discretisation
            # R = -(b[i] + s) #residual, function - boundary conditions + summed terms from discretisation
            # d = 0.5*omega*R #correction to current guess
            prevY = y[i]
            y[i] =  omega * Ynew/A[i,i]
            res[i] = np.abs(prevY - y[i]) #absolute value of residual, to check with tolerance

        iteration+=1
        if (iteration%500 == 0):
            print(iteration, np.max(res*y))
        else:
            pass
 return y
"""

#@jit(float64[:](float64[:,:],float64[:], float64[:], float64, float64), nopython = True)
@jit(nopython = True)
def new(A,x,f, bc1, bc2):
    rel_tol = 1e-6
    iteration = 0 
    omega = 1.00
    r0 = 1e-3
    h = dr
    #omega = 2.0 / (1.0 + np.sin(np.pi * h))
    l = len(x)
    #h = 1.0/(l-1)

    res = np.ones(len(x))
    b = np.zeros(len(x))
    b[0] = bc1
    b[-1] = bc2
    #tick = time.time()
    while (np.any(res > np.abs(np.multiply(rel_tol,x)))):
        for i in range(0, l):
            s = np.dot(A[i,:], x[:])
            xnew = f[i]*h**2 - (s - A[i,i]*x[i]) - b[i]
            xold = x[i]
            x[i] = xnew*omega/A[i,i]
            res[i] = np.abs(xold - x[i])
        iteration +=1
        if (iteration%10000 ==0):
            #print("\r"+ str(iteration),end='')
            print(iteration)
        else:
            pass 
    #print("\r"+ str(iteration))
    print(iteration)
    inv = np.linalg.inv(A)
    ans = np.dot(inv,f*h*h-b)
    #tock = time.time()
    #print(tick - tock)
    return x, ans 



"""sol = new(A,y,f, 0.0, 0.0)
#sol2 = th(A,y,f)
s = np.linalg.inv(A)
x = np.dot(s,f)
print(sol)
print(x)"""

@jit(nopython = True)
def new1(A,x,f, bc1, bc2,r):
    rel_tol = 1e-6
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
    b[0] = bc1
    b[-1] = bc2
    while (np.any(res > np.abs(np.multiply(rel_tol,x)))):
        for i in range(0, l):
            s = np.dot(A[i,:], x[:])
            xnew = -f[i]*h**2  - (s - A[i,i]*x[i])  - b[i]
            xold = x[i]
            x[i] = xnew*omega/A[i,i]
            res[i] = np.abs(xold - x[i])
        iteration += 1
        if (iteration%1000 == 0):
            #print("\r"+ str(iteration),end='')
            print(iteration)
        else:
            pass 
    #print("\r"+ str(iteration))
    print(iteration)
    """Determining answer via matrix inversion for comparison"""
    inv = np.linalg.inv(A)
    ans = np.dot(inv,-f*h*h - b)
    return x, ans 


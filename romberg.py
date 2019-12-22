import numpy as np 
import matplotlib.pyplot as plt

def romberg(f,a,b,err,iter_max):
    R = np.zeros((iter_max,iter_max), float)
    if (f(a) == 0.0) and (f(b) == 0.0):
        i = 1
    else:
        i = 0
    def trap_rule(f,a,b,N):
        h = (b - a)/N
        x = np.linspace(a,b,N+1)
        f_arr = f(x)
        sum_int = 0.0
        sum_int = np.sum(f_arr[1:-1])
        trap = 0.5*h*(f_arr[0] + f_arr[-1]) + h*sum_int
        return trap
    while i < iter_max :
        print(i)
        N = 2**i
        R[i,0] = trap_rule(f,a,b,N)
        for j in range(1,i+1):
            R[i,j] = (1.0/(4**j -1))*(R[i,j-1]*4**j - R[i-1,j-1])
        if (i>0):
            if(abs(R[i,j] - R[i,j-1]) < err ):
                break
        i +=1
    return R[i,j], i 
"""
def func(x):
    if (x > 1.0):
        y = 0.0
    else:
        y= x**4
    return y

def func2(x):
    return np.exp(-x*x)

A1 = romberg(func,0.0,2.0,1e-10,10)
print(A1)"""

def romberg_samp(y,x):
    length = len(y)
    a = x[0]
    b = x[-1]
    check = np.copy(length) - 1
    i = 0
    while check != 1 :
        check = check /2
        i += 1
    R = np.zeros((i,i), float)
    def trap_rule(y,x):
        N = len(x) -1
        h = (b - a)/N
        sum_int = 0.0
        sum_int = np.sum(y[1:-1])
        trap = 0.5*h*(y[0] + y[-1]) + h*sum_int
        return trap
    def find_nearest(arr,val):
        ind = (np.abs(arr[:] - val)).argmin()
        return ind
    for j in range(0, i):
        x_iter = np.linspace(x[0], x[-1], num = 2**j + 1, endpoint = 'true')
        y_iter = np.zeros(len(x_iter))
        for p in range(0, len(x_iter)):
            ind = find_nearest(x,x_iter[p])
            y_iter[p] = y[ind]       
        R[j,0] = trap_rule(y_iter,x_iter)
        for k in range(1,j+1):
            R[j,k] = (1.0/(4**k -1))*(R[j,k-1]*4**k - R[i-1,k-1])
    return R[j,k]
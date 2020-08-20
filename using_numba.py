import numpy as np 
from numba import jit , njit


@njit
def func(x, z,b,a):   
    if b ==5:
        print('number is right')
    out1 = a*x
    out2 = b*z
    return out1, out2
    
x = np.arange(1.0,5.0,step = 1.0)
z = np.ones(10)
a = np.ones(len(x))*2.0
b = 5
y1,y2 = func(x,z,b,a)
print(y1,y2)
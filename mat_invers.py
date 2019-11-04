import numpy as np
from numpy.linalg import inv
import matplotlib.pyplot as plt

def m_inv_lap(f,bc1,bc2):
    l = len(f)
    h = 1.0/(l-1)
    b = np.zeros(l)
    b[0] = bc1 
    b[-1] = bc2
    A = np.zeros((l,l))
    for i in range(0,l):
        A[i,i] = 2.0
        f[i] *= h**2
    for i in range(0,l-1):
        A[i+1,i] = -1.0
        A[i,i+1] = -1.0
    A_inv = inv(A)
    RHS = f + b
    y = A_inv @ RHS
    return y
"""f = np.zeros(4)
bc1 = 0.0
bc2 = 5.0

ysol = m_inv_lap(f, bc1, bc2)
ysol = np.insert(ysol, 0,bc1)
ysol = np.append(ysol, bc2)
print(ysol)
plt.figure()
plt.plot(np.arange(0,6,1.0),ysol)
plt.show()"""

import numpy as np
import iterative_sol as SOR
import discret_mat
import matplotlib.pyplot as plt 

l = 500
x = np.linspace(0.9,6.1,l)
phi = np.zeros(l)
phi0 = 1000.0
A = discret_mat.discret(x)
n = np.ones(l) * 10**19
e = -1.6*10**(-19) 
eps0 = 8.854*10**(-12)
dens0 = 10**19
x0 = 10**(-3)

phi = np.zeros(l)
phi[0], phi[-1] = -1.0, 0.0

n_n = n/dens0 

K = (x0**2 * dens0*e)/(phi0*eps0)

phi = SOR.SOR(A,phi/K,n_n,x)

fig, ax  = plt.subplots()
ax.plot(x,phi*phi0)
plt.show()
import numpy as np
import iterative_sol as SOR
import discret_mat
import matplotlib.pyplot as plt

x = np.arange(0.0,10.1,0.1)
phi0 = 10.0
bc1 = phi0
bc2 = phi0 
dens = (10**12)*np.ones(len(x))
eps0 = 8.854*10**(-12)
e = 1.6*10**(-19)

phi = np.zeros(len(x))
phi[0] = bc1
phi[-1] = bc2 

A = discret_mat.discret(x)
f = -e*dens/eps0
phi = SOR.SOR(A,phi,f,x) 


fig, ax = plt.subplots()
ax2 = ax.twinx() 
ax.plot(x, phi, 'blue')
ax2.plot(x,dens, 'red')
plt.show()
#Now the normalised version 

x_n = x/x[-1]
dens_0 = 10**9
dens_n = np.ones(len(x))
phi_n = np.zeros(len(x))
phi_n[-1] = 1.0 
phi_n[0] = 1.0 
K = (x[-1]**2)*e*dens_0/(eps0*phi0)
#phi_n[:] /= K 
f_n = -np.ones(len(x))
phi_n = SOR.SOR(A,phi_n, K*f_n,x_n)



#phi_n[:] *= phi0

fig, ax = plt.subplots() 
ax2 = ax.twinx() 
ax.plot(x_n, phi0*phi_n, 'blue')
ax2.plot(x_n,dens_n, 'red')
plt.show()
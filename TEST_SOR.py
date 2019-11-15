import numpy as np 
import iterative_sol 
import matplotlib.pyplot as plt 
import discret_mat

"A few key parameters for the solver"
l = 100 # grid length
x = np.linspace(0,1,l) # the grid 

dens = x**2 + np.exp(-x)  # test charge density
phi = np.zeros(l) 
phi[0] = 2.0 #boundary condition 1
phi[-1] = -1.0 #boundary condition 2
A = discret_mat.discret(x) # discretisation matrix 

phi = iterative_sol.SOR(A,phi,-dens,x)

phi_an = (1.0/12.0)*x**4 + np.exp(-x) - x*(25.0/12.0 + np.exp(-1)) + 1.0
fig, ax  = plt.subplots() 
l1 = ax.plot(x,phi, label = r'$\phi$', color = 'firebrick')
ax2 = ax.twinx()
l2 = ax2.plot(x,dens, label = r'$n_e$', color = 'navy')
l3 = ax.plot(x,phi_an, label = r'$\phi_{\mathrm{real}}$', color = 'green', linestyle = '--')
ax.set_xlabel(r'$x$')
ax.set_ylabel(r'$\phi$', rotation = 0)
ax2.set_ylabel(r'$n_e$', rotation = 0)
lines = l1 + l2 +l3
labs = [l.get_label() for l in lines]
ax.legend(lines, labs, loc = 0)
plt.show()
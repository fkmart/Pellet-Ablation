import numpy as np 
import matplotlib.pyplot as plt
import iterative_sol
import discret_mat

savedir = '/home/kyle/Documents/thesis/Local_figures'
title = 'polySORtest'
"""Integration constants"""


l = 100
x = np.linspace(0.0,1.0, l)

dens = 3.0*x**2 - 2.0*x

phi = np.zeros(l)
phi[0] = 0.0
phi[-1] = 0.0
A = discret_mat.discret(x)
phi = iterative_sol.SOR(A,phi, dens, x)

#a test on the solver is to differentiate the result 

check_the_first = np.dot(A,phi)
check_the_first /= (1.0/l)**2 # tested on 26/08/2019 and it works!!!

fig, ax = plt.subplots()
l1 = ax.plot(x, phi, color = 'royalblue', label = r'$\tilde{\phi}$')
ax2 = ax.twinx()
l2 = ax2.plot(x, dens,color = 'firebrick', label = r'$\tilde{\rho}_c = 3\tilde{x}^2 - 2\tilde{x}$')
ax.set_ylabel(r'$\tilde{\phi}$', fontsize = 12, rotation = 0)
ax.set_xlabel(r'$\tilde{x}$', fontsize = 12, rotation = 0)
ax2.set_ylabel(r'$\tilde{\rho}_c$', fontsize = 12, rotation = 0)

ax.xaxis.set_label_coords(0.50, -0.05)
ax.yaxis.set_label_coords(-0.04, 0.48)
ax2.yaxis.set_label_coords(1.04, 0.51)

lns = l1 + l2
labels = [l.get_label() for l in lns]

plt.legend(lns, labels, ncol = 1, loc = 'upper left')
plt.savefig(savedir+'/'+title+'.png', format = 'png', dpi = 1200)
plt.show()
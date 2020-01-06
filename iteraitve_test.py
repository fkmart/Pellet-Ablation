import numpy as np 
import matplotlib.pyplot as plt
import iterative_sol
import discret_mat
import os

#savedir = '/home/kyle/Documents/thesis/Local_figures'
savedir = os.path.join(os.getcwd(), 'Local_figures') + os.sep
title = 'SOR_poly'
"""Integration constants"""


l = 100
x = np.linspace(0.0,1.0, l)

#dens = np.sin(3.0*np.pi*x)
#dens = np.ones(l)
dens = -3.0*x**2 + 0.5*x - 1.0

phi = np.zeros(l)
phi[0] = 2.0
phi[-1] = 3.0
A = discret_mat.discret(x)
phi = iterative_sol.SOR(A,phi, dens, x)

#analytical check
#an = -(1.0/(9.0*np.pi**2))*np.sin(3.0*np.pi*x) - x
#an = 0.5*x**2 - 0.5*x + 1.0
#an = 0.5*x**2 + 2.5*x + 2.0
an = -0.25*x**4 + (1.0/12.0)*x**3 - 0.5*x**2 + (10.0/6.0)*x + 2.0
#a test on the solver is to differentiate the result 

check_the_first = np.dot(A,phi)
check_the_first /= (1.0/l)**2 # tested on 26/08/2019 and it works!!!

fig, ax = plt.subplots()
l1 = ax.plot(x, phi, color = 'royalblue', label = 'SOR Solution')
l3 = ax.plot(x,an, color = 'limegreen', linestyle = '--', label = r'$\tilde{\phi} = -\frac{\tilde{x}^4}{4} + \frac{\tilde{x}^3}{12} - \frac{\tilde{x}^2}{2} + \frac{5\tilde{x}}{3} + 2$')
ax2 = ax.twinx()
l2 = ax2.plot(x, dens,color = 'firebrick', label = r'$\tilde{n}_c = -3\tilde{x}^2 + \frac{1}{2}\tilde{x} -1$')
ax.set_ylabel(r'$\tilde{\phi}$', fontsize = 12, rotation = 0)
ax.set_xlabel(r'$\tilde{x}$', fontsize = 12, rotation = 0)
ax2.set_ylabel(r'$\tilde{n}_c$', fontsize = 12, rotation = 0)

ax.xaxis.set_label_coords(0.50, -0.05)
ax.yaxis.set_label_coords(-0.04, 0.43)
ax2.yaxis.set_label_coords(1.04, 0.49)

lns = l1 + l2 + l3
labels = [l.get_label() for l in lns]

plt.legend(lns, labels, ncol = 1, loc = 'lower center')
plt.savefig(savedir + title+'.png', format = 'png', dpi = 1200)
plt.show()
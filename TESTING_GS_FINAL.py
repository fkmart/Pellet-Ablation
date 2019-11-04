import numpy as np
import discret_mat
import iterative_sol
import matplotlib.pyplot as plt
import scipy.integrate as spinteg

from TEST import new 
#from phi_calc_total import phi_tot
import charge_dens
import elec_phi_calc 

#constants for electrostatics
eps = 8.854*10**(-12)
e = 1.6*10**(-19)

con = 10**10
K = (e/eps)*con


"""Variables for the integration"""
l = 101
x = np.linspace(0.0,1.0,l)
rhoc_e = np.zeros(l) + K
rhoc_i = np.zeros(l)
e_dens = 5.0
off = 0
#rhoc_e[int(l/2)-1 -off : int(l/2)+1 - off] = -e_dens

tot_e = spinteg.simps(rhoc_e[int(l/2)-1 -off : int(l/2)+1 - off], x[int(l/2)-1 -off : int(l/2)+1 - off])

fac = 1.0
thick = fac*np.count_nonzero(rhoc_e)
rhoc_i[-int(thick)+int(fac-1):l] = e_dens/fac

tot_i = spinteg.simps(rhoc_i[-int(thick)+int(fac-1):l], x[-int(thick)+int(fac-1):l])
"""Checking the charge numbers match"""
check = tot_e/tot_i
print(check)
#if((np.abs(np.abs(check) - 1.0)) > 1e-7):
#    raise ValueError("Ion and electron numbers don't balance")

bc1 = 0.0
bc2 = 0.0

"""Electrons"""
phi_e = np.zeros(l)
phi_e[0] = bc1
phi_e[-1] = bc2
A = discret_mat.discret(x)
phi_e, phi_mat_e = iterative_sol.SOR(A, phi_e, rhoc_e, x)


"""Ions"""
phi_i = np.zeros(l)
phi_i[0] = bc1
phi_i[-1] = np.amax(np.abs(phi_e))
#phi_i[-1] = bc2
phi_i, phi_mat_i = iterative_sol.SOR(A, phi_i, rhoc_i, x)

"""Totals"""

phi_t = np.zeros(l)
phi_t[0] = bc1
phi_t[-1] = phi_i[-1]
rhoc_t = rhoc_e + rhoc_i
phi_t, phi_mat_t = iterative_sol.SOR(A, phi_t, rhoc_t, x)

"""Plotting"""

fig = plt.figure()
plt.grid()
ax = fig.add_subplot(1,1,1)
ax2 = ax.twinx()
l1 = ax.plot(x,phi_e, label = r'$\phi_e$')
#l2 = ax.plot(x, phi_i*K, label = r'$\phi_i$')
#l3 = ax.plot(x, phi_i + phi_e, label = r'$\phi_e + \phi_i$')
l4 = ax2.plot(x,rhoc_e, 'k--',label = r'$\rho_e$')
#l5 = ax2.plot(x, rhoc_i*con, 'r-.'  , label = r'$\rho_i$')
#l6 = ax.plot(x, phi_t*K, label = r'$\phi_t$')

lns = l1# + l2 + l4+ l5 + l6
labels = [l.get_label() for l in lns]


"""Analytical result to compare"""
#phi_test = -0.5*rhoc*x**2 +(0.5*rhoc*x[-1] + phi[-1]/x[-1])*x + phi[0]

#plt.plot(x,phi_test, '--', label = 'Analytical')
#plt.plot(x,phi_mat)

"""Testing the 'point' charge approach by spreading the charge over a thin layer for each species
and testing how the potential looks. """
ax.set_xlabel(r' $\tilde{r}$')
ax.set_ylabel(r'$\tilde{\phi}$', rotation = 0)
ax2.set_ylabel(r'$\tilde{\rho_c}$', rotation = 0)
plt.legend(lns, labels, ncol = 2, loc = 'upper left')
plt.show()
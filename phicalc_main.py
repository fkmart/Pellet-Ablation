
from gen_var import dr, rc, rp, r0, e
import os
import scipy.integrate as spinteg
import scipy.interpolate as spint
import number
from MB_calc import MB
import numpy as np
import matplotlib.pyplot as plt
import charge_dens
import elec_phi_calc
import discret_mat
import iterative_sol

t = 2 #  remove this line from code when incorporated into main code via "import" statement

e_dens, r, norm, tot = charge_dens.e_dens(t)

N = number.number(t)

"""SEttign up the ion density given total number of particles recently injected into the system"""
i_dens = np.zeros(len(r))
i_dens[-1] = N

"""boundary conditions"""
phi_p = 0.0
phi_c = 0.0
phi = np.zeros(len(r))

phi[0] = phi_p
phi[-1] = phi_c

A = discret_mat.discret(r)

i_dens = np.zeros(len(r))
i_dens[-1] = tot

q_dens = e_dens + i_dens
phi = iterative_sol.SOR(A, phi, q_dens, r)

#phi, phi_mat = elec_phi_calc.elec_phi_test(r, q_dens)

plt.plot(r,phi)
plt.show()
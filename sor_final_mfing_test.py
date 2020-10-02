# -*- coding: utf-8 -*-
"""
Created on Thu Sep 24 09:16:20 2020

@author: Kyle
"""
import numpy as np 
import iterative_sol as SOR 
import matplotlib.pyplot as plt 
import discret_mat as dm

n0 = 10**19 
eps0 = 8.854*10**(-11)
x0 = 1e-3
e = 1.6*10**(-19)
phi_0 = 100.0

l = 257
x_n = np.linspace(-1.0, 1.0, num = l)
n_n = np.ones(l) 
phi_n = np.zeros(l)

phi_n[0] = -1000.0
phi_n[0] /= phi_0

A = dm.discret(x_n)
phi_out = SOR.SOR(A,phi_n*phi_0,e*n0*n_n/eps0,x0*x_n)

K = e*x0*x0*n0/(eps0*phi_0)


phi_n = SOR.SOR(A, phi_n, K*n_n, x_n)

fig, ax = plt.subplots()
ax.plot(x_n,phi_out)
ax.plot(x_n,phi_n*phi_0)
plt.show()
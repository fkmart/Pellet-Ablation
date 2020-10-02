# -*- coding: utf-8 -*-
"""
Created on Fri Aug 21 16:07:22 2020

@author: Kyle
"""
import numpy as np
import matplotlib.pyplot as plt 

import romberg as ro
def solver(n_old, D,h,delta_t, n_new):
    for j in range(1, len(n_old) - 1):
        diff1 = n_old[j+1] - n_old[j]
        diff2 = n_old[j] - n_old[j-1]
        term1 = D[j] * (diff1)/h
        term2 = D[j-1] * (diff2)/h
        n_new[j] = n_old[j] + (delta_t/h)*(term1 + term2)
    return n_new 


def solver2(n_old, D, h, delta_t, n_new):
    for i in range(1, len(n_old) - 1):
        term1 = D[i] *(n_old[i+1] + n_old[i-1] - 2.0*n_old[i])/(h*h)
        term2 = ((D[i+1] - D[i])/h)*(n_old[i+1] - n_old[i])/h
        n_new[i] = n_old[i] + delta_t *(term1 + term2)
    return n_new
l = 129

x = np.linspace(-1.0,1.0, num = l)
il = 240
gap = 50
iu = il + gap
dens = np.zeros(l)
dens = np.linspace(0.0, 2.0, num = l)
#dens = np.sin((np.pi*x/6.0))
dens[0] = 0.0
dens[-1] = 0.0
#dens[il:iu] = 1.0

check = ro.romberg_samp(dens,x)
print('Original mass is ' + str(check))
D = np.linspace(0.01, 0.1, num = l)*0.005
#D = np.ones(l)*0.0005
#D[:] = 1e-2
dt = 0.1
h = x[1] - x[0]

print('r is ' + str(D[0]*dt/(h*h)))

n_new = np.zeros(l)
n_new[0] = 0.0
n_new[-1] = 0.0
n_out = np.copy(dens)

fig,ax = plt.subplots()
ax.plot(x, dens, label = r'$\Delta t = 0$')
plt.xlabel(r'$\tilde{r}$', fontsize = 12, rotation = 0)
plt.ylabel(r'$\tilde{n}$', fontsize = 12, rotation = 0)
plt.yscale('log')
for a in range(1, 100):
    w = a-1
    n_out = solver2(n_out, D, h, dt, n_new)
    checker = ro.romberg_samp(n_out,x)
    print('Mass after diffusion is ' + str(checker))
    if a%20 ==0:
        ax.plot(x, n_out, label = r'$\Delta t = $' + str(a))
plt.legend()
#plt.savefig('FD_Diffusion.png', format = 'png', dpi = 1200)
plt.show()

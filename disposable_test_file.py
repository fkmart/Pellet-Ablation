# -*- coding: utf-8 -*-
"""
Created on Tue Aug 18 19:22:25 2020

@author: Kyle
"""

import numpy as np
import matplotlib.pyplot as plt 
import romberg as ro
import scipy.special as scisp


l = 1025

x = np.linspace(-5.0,5.0, num = l)
dens = np.zeros(l)

mid = int(l*0.5)

dens[300] = 1.0

cs = 1e-20
nd = 1e23
mfp = (cs*nd)**-1
vth = 1e5


D = mfp*vth
t = 1e-5
B = ro.romberg_samp(dens,x)
diff_dens= (B/(np.sqrt(4.0*np.pi*D*t)))*np.exp(-(x - x[300])**2 /(4*D*t))
#diff_dens = 0.5*(dens[mid]/(2.0*np.sqrt(np.pi*D*t)))*scisp.erfc(x**2/(4*D*t))
C = ro.romberg_samp(diff_dens,x)
print(B)
print(C)
fig,ax = plt.subplots()
ax.plot(x,dens)
ax.plot(x,diff_dens)

plt.show()
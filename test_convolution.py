# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 13:19:01 2020

@author: Kyle
"""

import numpy as np 
import scipy.signal as scisig 
import matplotlib.pyplot as plt
import romberg as ro
import os 

direc = os.getcwd()

load_dir = os.path.join(direc,'many_iteration','neutral','500eV','analysed_outputs') + os.sep 

file = np.genfromtxt(load_dir + 'outputs1.txt', delimiter = ',', dtype = 'str')

dens = file[5:,1]
dens = np.asarray([float(w) for w in dens])
r = file[5:,0]
r = np.asarray([float(w) for w in r])
r_cent = (r[0] + r[-1]) * 0.5

rl = 25
rr = 726

#dens = dens[rl:rr]
#r = r[rl:rr]

r_cent = (r[-1] + r[0])*0.5

#get diffusion properties 
neut_dens = 0.01*(1.0 + r[25]**2)/(1.0 + r**2)
D = neut_dens**-1
#D = 50.0
t = 0.0005
prefac = 1.0/(np.sqrt(4.0*np.pi * D * t))
Z = (r - r_cent)/np.sqrt((4.0*D*t))
win = prefac*np.exp(-Z**2)

output = scisig.convolve(dens, win, mode = 'same')
A = np.sum(dens)
B = np.sum(output)
output *=A/B
output[:26] = 0.0 
output[701:] = 0.0
fig, (ax1,ax2,ax3) = plt.subplots(3,1, sharex = True)
ax1.plot(r,dens)
ax2.plot(r,win)
ax3.plot(r, output)
plt.show()

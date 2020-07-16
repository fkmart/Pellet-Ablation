# -*- coding: utf-8 -*-
"""
Created on Tue Jul  7 09:39:39 2020

@author: Kyle
"""

import numpy as np 
import matplotlib.pyplot as plt 
import stop_calc_rp_rc as scr
import os 

direc = os.getcwd() 
TYPE = ['neutral','poisson','sheath']
style = ['-',':','--']
col = ['blue','red','green']
energy = ['500eV','1000eV','1500eV']
t = np.linspace(0,19,num = 20)

fig, ax = plt.subplots()
types = TYPE[0]

load_dir = os.path.join(direc,'many_iteration','neutral', '1500eV','analysed_outputs') + os.sep
    #t_params = np.genfromtxt(load_dir + 't_params.txt',dtype = 'str')
    #t_numbs = [float(s) for s in t_params[1,:]]
    #t_arr = np.arange(t_numbs[0], t_numbs[1], t_numbs[2])
    #indices = np.loadtxt(load_dir + 'time_indices.txt')
rp_rc_arr = np.loadtxt(load_dir + 'rp_rc.txt')
rps = rp_rc_arr[:,0]
rcs = rp_rc_arr[:,1]
l = len(rps)
plt.xlim(rps[-1],4)
for k in range(0,l):
    rp = rps[k]
    rc = rcs[k]
    r = np.linspace(rp, rc, num = 1000, endpoint = 'true')
    dens = 0.01*(1.0 + rp**2)/(1.0 + r**2)
    ax.plot(r, dens)
    ax.set_yscale('log')
plt.show()
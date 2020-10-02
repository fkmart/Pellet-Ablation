# -*- coding: utf-8 -*-
"""
Created on Tue Sep  1 12:25:51 2020

@author: Kyle
"""
import numpy as np
import matplotlib.pyplot as plt 
import os 
from gen_var import many_start, rp_hr , t_hr, delta_t

direc = os.getcwd() 
TYPE = ['neutral','poisson','sheath'] 

l = 2000
rp_model = np.zeros((len(TYPE),l-1))

fig, ax = plt.subplots()
ax2 = ax.twinx()
rpd_ngs = rp_hr[many_start:]**(-2.0/3.0)
lines = []
lines = ax2.plot(rp_hr[many_start:], rpd_ngs,label = 'NGS model', color = 'black')
for i in range(0, len(TYPE)):
    load_dir = os.path.join(direc, 'many_iteration', TYPE[i], '1000eV', 'analysed_outputs') + os.sep
    for j in range(1, l, 50):
        file = np.genfromtxt(load_dir + 'outputs' + str(j) + '.txt', dtype = 'str', delimiter = ',')
        rp_model[i,j-1] = np.float(file[1,5])
    rp_uni = np.unique(rp_model[i,:])
    times = np.ones(len(rp_uni))
    delta_rp = np.zeros(len(rp_uni) - 1)
    for j in range(0, len(rp_uni) -1 ):
        times[j] = list(rp_model[i,:]).count(rp_uni[j])
        delta_rp[j] = rp_uni[j+1] - rp_uni[j]
    rpd = delta_rp[:]/(times[:-1]*delta_t)
    lines += ax.plot(rp_uni[1:-1],rpd[1:], label = TYPE[i])

labs = [d.get_label() for d in lines]
plt.legend(lines, labs)
plt.show()
        
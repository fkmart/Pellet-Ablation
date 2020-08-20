# -*- coding: utf-8 -*-
"""
Created on Sun Aug  9 08:10:19 2020

@author: Kyle
"""

import numba_test_rkf as rkf
import step_calc as sc
import numpy as np
from numba import njit 
from gen_var import t_hr, I, trunc_fac
from gen_var import rp_hr as rp
from gen_var import rc_hr as rc
from electron import RME, M_fac
from stopblock import stopblock_phi_mod_rkf_jit as func

BT = np.zeros((6,6))
BT[1,:] = [0.25,0.25,0.0, 0.0, 0.0, 0.0]
BT[2,:] = [0.125, 3.0/32.0, 9.0/32.0, 0.0, 0.0, 0.0]
BT[3,:] = [12.0/13.0, 1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, 0.0, 0.0]
BT[4,:] = [1.0, 439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0, 0.0]
BT[5,:] = [0.5, -8.0/27.0, 2.0, -3544.0/2565.0, 1859./4140.0 , -11.0/40.0]

#butcher-tableau for solutions

BTS = np.zeros((2,6))
BTS[0,:] = [25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, - 0.2, 0.0]
BTS[1,:] = [16.0/135.0, 0.0, 6656.0/12825.0,28561.0/56430.0, -9.0/50.0, 2.0/55.0 ]

I_non = I/(RME*M_fac)

i = 150000


E = np.arange(250,15000, step = 500)
@njit
def function(E,r,i,phi):
    energy_outputs = []
    space_outputs = []
    lengths = np.zeros(len(E))
    for j in range(0, len(E)):
        print(j)
        en = [E[j]/(RME*M_fac)]
        #k1,k2,k3,k4,k5,k6 = 0.0,0.0,0.0,0.0,0.0,0.0 #RKF4 algorithm terms
        #z = np.arange(0,lr-1,1) #helps to determine r index in while loop with y counter
        #y = 0 #reset y counter to 0 after every energy
        xr = rc[i]
        xl = rp[i]
        dphi = 0.0
        tot_pot = 0.0
        x = rc[i]
        x_tot = [x]
        I_non = I/(RME*M_fac)
        y_steps = np.zeros(5)
        dr = 1e-3
        c = 0
        while en[c] > trunc_fac*I_non and en[c] - y_steps[0] >trunc_fac*I_non and en[c] -y_steps[1] > trunc_fac*I_non and en[c] - y_steps[2] > trunc_fac*I_non and en[c] - y_steps[3] > trunc_fac*I_non and en[c] - y_steps[4] > trunc_fac*I_non and x > rp[i] and en[c] + dphi > trunc_fac*I_non:
            x,energy,h,err, k = rkf.rkf_jit(x,en[-1],dr,BT,BTS,xr,xl,I_non)
            #k1,k2,k3,k4,k5,k6 = new_k[0], new_k[1], new_k[2], new_k[3], new_k[4], new_k[5]
            en.append(energy)
            c +=1
            x_tot.append(x)
            for p in range(0,5):
                y_steps[p] = sc.ys_calc_jit(h,k,BT,p)
        else:
            en = en[:] 
            x_tot = x_tot[:]
            tot_pot *= RME*M_fac
            for s in range(0, len(en)):
                en[s] *= RME*M_fac
                energy_outputs.append(en[s])
                space_outputs.append(x_tot[s])
            lengths[j] = len(en)
            #ener_prof = np.asarray(en)
            #ener_prof = en.append(x_tot) # add the spacial points to ener_prof as final output
            #ener_prof = np.reshape(ener_prof, (int(len(ener_prof)/2),2), 'F') # reshape array for two columns                  
            #output_list = output_list.append(ener_prof)
    print('done')
    return energy_outputs,space_outputs, lengths
    
r = np.linspace(rp[i], rc[i], num = 10000)
phi = np.zeros(10000)
start = 0
out1, out2, lengths = function(E,r,i,phi)
for h in range(0,len(lengths)):
    reform = np.asarray(out1[start:int(lengths[h])])
    reform2 = np.asarray(out2[start:int(lengths[h])])
    saving = np.append(reform,reform2)
    saving = np.reshape(saving, (int(lengths[h]),2), 'F')
   
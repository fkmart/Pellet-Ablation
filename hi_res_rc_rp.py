# -*- coding: utf-8 -*-
"""
Created on Fri Jul  3 09:51:55 2020

@author: Kyle
"""

import stop_calc_rp_rc as sc 
#from gen_var import dt_hr, t_low_hr, t_upper_hr, t_hr 
import numpy as np 
import math as mt
import rc_solve


dt_hr = 1e-6
t_low_hr = 1e-4
t_upper_hr = 0.5
t_hr = np.arange(t_low_hr,t_upper_hr,dt_hr, endpoint = 'true')
rp_hr = sc.calc_rp(t_hr) 
np.save('rp_hr.npy', rp_hr)
rc_hr = rc_solve.rc_sol(t_hr) 
np.save('rc_hr.npy', rc_hr)
rpd_hr = np.zeros(len(t_hr)) 

for i in range(0, len(t_hr)):
    rpd_hr[i] = -0.5*(1.0/(np.sqrt(np.log(mt.cosh(1)))))*(mt.tanh(1 - t_hr[i])/(np.sqrt(np.log(mt.cosh(1 - t_hr[i])))))
    
np.save('rpd_hr.npy', rpd_hr)
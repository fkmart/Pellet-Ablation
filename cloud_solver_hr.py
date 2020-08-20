# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 15:58:57 2020

@author: Kyle
"""

import numpy as np 

t = np.arange(0.0,0.5, step = 1e-6)

import stop_calc_rp_rc as scr 

rp_hr = scr.calc_rp(t)
rc_hr = scr.calc_rc(t)

np.savetxt('rp_hr.txt', rp_hr)
np.savetxt('rc_hr.txt', rc_hr)
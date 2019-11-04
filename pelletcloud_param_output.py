#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 12:50:16 2019

@author: kyle
"""

from stop_calc_rp_rc import *
import stop_calc_rp_rc 
import numpy as np

lt = 10
t = np.linspace(0.09, 0.99, 90)
rp = stop_calc_rp_rc.calc_rp(t)
rc = stop_calc_rp_rc.calc_rc(t)
rcdot = stop_calc_rp_rc.rcdot(t, rc, rp)
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan  4 15:08:37 2019

@author: kyle
"""

import numpy as np
from gen_var import *


def stop(x,E):
    #den = eps*dens*(1.0 +(rp[j*int(len(rp)/lt)])**2)/(1.0 + (rc[j*int(len(rc)/lt)]-x)**2)
    nondim = r0/RME
    c1= 0.153536
    dxde_gas = (nondim*den*(np.log((E**2.0)*(E + 2.0)/2.0) + (1.0 + (E**2.0)/8.0 -
               (2.0*E +1.0)*np.log(2.0))/(E + 1.0)**2.0 - 2.0*np.log(I/(RME*M_fac)) 
               - delta)*zovera*(-1.0*(c1/(2.0*E))))
    return dxde_gas 


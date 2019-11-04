#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 11:48:02 2019

@author: kyle
"""

import numpy as np

def sph_integ_int(integrand, r,i):
    integ = integrand[0:i]
    integ = integ*r0**(3) # for integrating over density
    out = []
    for j in range(0, i):
        y = (((4*np.pi)/3.0)*(r[j]**3 - r[j-1]**3)*integ[j])/(4.0*np.pi*r[j]**2)
        out = np.append(out,y)
    out_sum = np.sum(out)
    return out_sum

def sph_integ_ext(integrand, r, i):
    #integ = integrand[i, len(r)]
    integrand = integrand*r0**3 # for integrating over density
    out = []
    for j in range(i, len(r)):
       y = ((4.0/3.0)*np.pi*(r[j]**3 - r[j-1]**3)*integrand[j])/(4.0*np.pi*(r[j]- r[i])**2)    
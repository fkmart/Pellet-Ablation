# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 11:25:36 2020

@author: Kyle
"""

import numpy as np 
import rkf45_sheath as rkf 

def sheath_solver():
    from electron import e_bar as phi_0
    """Real units"""
    e = 1.6*10**(-19) # in C
    kb = 1.38*10**(-23) # in J/K
    mi = 1.67*10**(-27) #in kg
    n0 = 10**19 #in m^-3
    eps0 = 8.854*10**(-12) # in F/m
    me = 9.11*10**(-31) # in kg
    #phi_0 = 10**3 # in v
    
    
    T = phi_0*e/kb # in K
    u0 = np.sqrt(2*kb*T/mi) # in m/s
    ld = (eps0*kb*T/(e*e*n0))**(0.5) # in m
    epsp = mi*u0*u0/(2.0*e) #in V
    
    #FInal constant
    K = (2.0*n0*e*ld*ld/(phi_0*eps0))**(-0.5)
    
    """Normalised Units"""
    
    T_n = T*kb/e #in V
    T_n /= phi_0 #non-dim
    epsp_n = epsp/phi_0 #non-dim
    e_n = e/e #non-dim
    n_n = n0/n0 # non-dim
    """Setting up the integrand matrix"""
    phi_s = np.log(np.sqrt(2.0*mi/(2.0*np.pi*me)))
    #phi_s /=phi_0
    l = 3000
    phi_1 = -1e-3
    #phi = np.linspace(phi_1, phi_s, l)
    #phi = np.logspace(np.log10(phi_1), np.log10(phi_s), l)
    x0 = 0.0
    h = 0.01
    
    #epsp_n *=phi_s*1.01
    
    phi, x, dr_arr, err_arr = [],[],[],[]
    phi = [-phi_s]
    x = [x0]
    """Setting up the function to be integrated"""
    count = 0
    def sheath_final(phi_n,x):
        term1 = T_n*np.exp(phi_n/T_n) - T_n
        term2 = 2.0*epsp_n*(1.0 - phi_n/epsp_n)**(0.5) - 2.0*epsp_n
        diff1 = term1 + term2 
        diff = K*(diff1)**(-0.5)
        return diff
    """Normalised integration"""
    while phi[-1] < phi_1:
        phi_out, x_out, dr_out, err_out = rkf.rkf(phi[-1],x[-1],h,sheath_final,-phi_1,-phi_s)
        phi = np.append(phi, phi_out)
        x = np.append(x, x_out)
        dr_arr = np.append(dr_arr, dr_out)
        err_arr = np.append(err_arr, err_out)
        count +=1
    return x[-1]* ld, phi*phi_0,x*ld
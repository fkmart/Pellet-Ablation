# -*- coding: utf-8 -*-
"""
Created on Thu Jun 25 13:03:58 2020

@author: Kyle
"""

import numpy as np 
import rk45_mod as rk
#from gen_var import eps
def rc_sol(t):
    dt = t[1] - t[0]
    c = 1.0
    def func_rc(t,rc):
     eps = 0.01
     drcdt = ((eps - 1.0)/eps)*((((np.log(np.cosh(c*(1 - t)))/(np.log(np.cosh(c))))**(0.5))**2)/rc**2)*((1.0 
     # ((eps - 1.0)/eps)*((rnorm[i]**2)/(r**2))*((1.0 +
     + rc**2)/(1.0 +((np.log(np.cosh(c*(1 - t)))/(np.log(np.cosh(c))))**(0.5))**2))*(-0.5*c)*((np.log(np.cosh(c)))**(-0.5))*(((np.log(np.cosh(c*(1.0 
     # r**2)/(1.0 + rnorm[i]**2))*drdtnorm[i]
     - t)))))**(-0.5))*(np.tanh(c*(1-t))) - 2.0*(((np.log(np.cosh(c*(1.0 - t)))/(np.log(np.cosh(c))))**(0.5))*-0.5*c*((np.log(np.cosh(c)))**(-0.5))*(((np.log(np.cosh(c*(1.0 
     - t)))))**(-0.5))*(np.tanh(c*(1-t)))/rc**2)*((1.0 + rc**2)/(1.0 +(np.log(np.cosh(c*(1.0 - t)))/(np.log(np.cosh(c))))**(1.0)))*(rc - 
     (np.log(np.cosh(c*(1 - t)))/(np.log(np.cosh(c))))**(0.5) + np.arctan((((np.log(np.cosh(c*(1.0 
     - t)))/(np.log(np.cosh(c))))**(0.5) -rc)/(1.0 +rc*(np.log(np.cosh(c*(1 - t)))/(np.log(np.cosh(c))))**(0.5)))))
     return drcdt 
    
    rc_out = []
    time = t[0]
    rc = [1.0]
    for i in range(0,len(t)):
        rc_out = rk.rk(time,rc[-1], func_rc, dt)
        time += dt
        rc = np.append(rc, rc_out)
    return rc[1:]
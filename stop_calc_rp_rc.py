#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 27 10:45:28 2019

@author: kyle
"""

import numpy as np
import scipy.integrate as spit
import matplotlib.pyplot as plt

k = 1.0
k1 = 1.0
k2 = 1.0

c1 = 1.0
c2 = 1.0
c3 = 1.0

t = np.linspace(0, 0.999, 999)

rc0 = np.zeros(1)
rc0[0] = 1.0

eps = 0.01


#dt = t[1] - t[0]
dadt = []
drdtnorm =[]
rnorm = []
def calc_rp( t):
    lt = len(t)
    c = 1.0
    #for i in range(0, len(t)):
     #   dadt[i] = c1*np.tanh(c2*(1 - c3*t[i]))
      #  rnorm[i] = k2*(np.log(np.cosh(k1*(1 - t[i])))/(np.log(np.cosh(k1))))**(0.5)
       # drdtnorm[i] = -0.5*k2*k1*((np.log(np.cosh(k1)))**(-0.5))*(((np.log(np.cosh(k1*(1 - t[i])))))**(-0.5))*(np.tanh(k1*(1-t[i])))
    dadt = c1*np.tanh((c*(1-c3*t)))
    rnorm = k2*(np.log(np.cosh(c*(1 - t)))/(np.log(np.cosh(c))))**(0.5)
    drdtnorm = -0.5*k2*c*((np.log(np.cosh(c)))**(-0.5))*(((np.log(np.cosh(c*(1 - t)))))**(-0.5))*(np.tanh(c*(1-t)))
    #np.save('rp_t'+str(lt)+'.npy', rnorm)
    return rnorm
    
#Setting up differential equation through a function definition

#SOlving for rc via 

def calc_rc( time):
    c = 1.0
    lt = len(time)
    dt = time[1] - time[0]
    fac = 1
    t = np.linspace(time[0], time[len(time) -1], fac*len(time))
    dt = t[1] - t[0]
    def func_rc(t,rc):
     
     drcdt = ((eps - 1.0)/eps)*((((k2*np.log(np.cosh(c*(1 - t)))/(np.log(np.cosh(c))))**(0.5))**2)/rc**2)*((1.0 
     # ((eps - 1.0)/eps)*((rnorm[i]**2)/(r**2))*((1.0 +
     + rc**2)/(1.0 +((k2*np.log(np.cosh(c*(1 - t)))/(np.log(np.cosh(c))))**(0.5))**2))*(-0.5*k2*c)*((np.log(np.cosh(c)))**(-0.5))*(((np.log(np.cosh(c*(1.0 
     # r**2)/(1.0 + rnorm[i]**2))*drdtnorm[i]
     - t)))))**(-0.5))*(np.tanh(c*(1-t))) - 2.0*((k2*(np.log(np.cosh(c*(1.0 - t)))/(np.log(np.cosh(c))))**(0.5))*-0.5*k2*c*((np.log(np.cosh(c)))**(-0.5))*(((np.log(np.cosh(c*(1.0 
     - t)))))**(-0.5))*(np.tanh(c*(1-t)))/rc**2)*((1.0 + rc**2)/(1.0 +k2*(np.log(np.cosh(c*(1.0 - t)))/(np.log(np.cosh(c))))**(1.0)))*(rc - 
     (np.log(np.cosh(c*(1 - t)))/(np.log(np.cosh(c))))**(0.5) + np.arctan((((np.log(np.cosh(c*(1.0 
     - t)))/(np.log(np.cosh(c))))**(0.5) -rc)/(1.0 +rc*(np.log(np.cosh(c*(1 - t)))/(np.log(np.cosh(c))))**(0.5)))))
     return drcdt 

    rcloud = spit.solve_ivp(func_rc, (t[0], t[len(t) -1]), rc0, method = 'RK45', max_step = dt)
    rc = rcloud["y"]
    #print(rc)
    rc = np.ravel(rc)
    rc = rc[:].copy()
    rc = rc[::fac]
    #print(rc)
    #np.save('rc_t'+str(lt)+'.npy', rc)
    return rc

def rcdot(t,rc, rp): #loop over the time array as we want to flow speed over the whoel cloud at a singluar point in time for all times 
    lt = len(t)
    out = []
    r_pl = np.zeros((lt, 500))
    for i in range(0, len(t)):
        
        def func_rc(t, rc):
         
         drcdt = ((eps - 1.0)/eps)*((((k2*np.log(np.cosh(k1*(1 - t)))/(np.log(np.cosh(k1))))**(0.5))**2)/rc**2)*((1.0 
         
         + rc**2)/(1.0 +((k2*np.log(np.cosh(k1*(1 - t)))/(np.log(np.cosh(k1))))**(0.5))**2))*(-0.5*k2*k1)*((np.log(np.cosh(k1)))**(-0.5))*(((np.log(np.cosh(k1*(1.0 
         
         - t)))))**(-0.5))*(np.tanh(k1*(1-t))) - 2.0*((k2*(np.log(np.cosh(k1*(1.0 - t)))/(np.log(np.cosh(k1))))**(0.5))*-0.5*k2*k1*((np.log(np.cosh(k1)))**(-0.5))*(((np.log(np.cosh(k1*(1.0 
         - t)))))**(-0.5))*(np.tanh(k1*(1-t)))/rc**2)*((1.0 + rc**2)/(1.0 +k2*(np.log(np.cosh(k1*(1.0 - t)))/(np.log(np.cosh(k1))))**(1.0)))*(rc - 
         (np.log(np.cosh(k1*(1 - t)))/(np.log(np.cosh(k1))))**(0.5) + np.arctan((((np.log(np.cosh(k1*(1.0 
         - t)))/(np.log(np.cosh(k1))))**(0.5) -rc)/(1.0 +rc*(np.log(np.cosh(k1*(1 - t)))/(np.log(np.cosh(k1))))**(0.5)))))
         return drcdt
        rcdot = func_rc(t[i], rc[i])
        densrc = eps*(1.0 + rp[i]**2)/(1.0 + rc[i]**2)
        constant = densrc*rcdot
        r = np.linspace(rp[i], rc[i], 500)
        dens = eps*(1.0 + rp[i]**2)/(1.0 + r**2)
        v_flow = constant/dens
        out.append(v_flow)
        r_pl[i,:] = r
    out = np.reshape(out, (len(out),len(r)), 'F')
    return out, r_pl

def rdot(t,rc,rp):
    lt = len(t)
    r_pl = np.zeros((lt,500))
    rcd = np.zeros(lt)
    out = np.zeros((lt,500))
    eps = 0.01
    densrc = np.zeros(lt)
    for i in range(0, lt):
        def func_rc(t, rc):
         
         drcdt = ((eps - 1.0)/eps)*((((k2*np.log(np.cosh(k1*(1 - t)))/(np.log(np.cosh(k1))))**(0.5))**2)/rc**2)*((1.0 
         
         + rc**2)/(1.0 +((k2*np.log(np.cosh(k1*(1 - t)))/(np.log(np.cosh(k1))))**(0.5))**2))*(-0.5*k2*k1)*((np.log(np.cosh(k1)))**(-0.5))*(((np.log(np.cosh(k1*(1.0 
         
         - t)))))**(-0.5))*(np.tanh(k1*(1-t))) - 2.0*((k2*(np.log(np.cosh(k1*(1.0 - t)))/(np.log(np.cosh(k1))))**(0.5))*-0.5*k2*k1*((np.log(np.cosh(k1)))**(-0.5))*(((np.log(np.cosh(k1*(1.0 
         - t)))))**(-0.5))*(np.tanh(k1*(1-t)))/rc**2)*((1.0 + rc**2)/(1.0 +k2*(np.log(np.cosh(k1*(1.0 - t)))/(np.log(np.cosh(k1))))**(1.0)))*(rc - 
         (np.log(np.cosh(k1*(1 - t)))/(np.log(np.cosh(k1))))**(0.5) + np.arctan((((np.log(np.cosh(k1*(1.0 
         - t)))/(np.log(np.cosh(k1))))**(0.5) -rc)/(1.0 +rc*(np.log(np.cosh(k1*(1 - t)))/(np.log(np.cosh(k1))))**(0.5)))))
         return drcdt
        rcd[i] = func_rc(t[i], rc[i])
    """ densrc[i] = eps*((1.0 + rp[i]**2)/(1.0 + rc[i]**2))
    for i in range(0, lt-1):
        av = (rcd[i]*densrc[i] + rcd[i+1]*densrc[i+1])*0.5
        r = np.linspace(rp[i], rc[i], 500)
        dens = eps*(1.0 + rp[i]**2)/(1.0 + r[:]**2)
        out[i,:] = av/dens[:]
        r_pl[i,:] = r[:]
    print(rcd[990])"""
    return rcd # out, r_pl 

def u_solution(t1,t2,l):
    from gen_var import rc,rp, eps
    from gen_var import t as time
    import math as mt 
    shift = t2-t1 
    rp1,rp2 = rp[t1], rp[t2]
    rc1,rc2 = rc[t1], rc[t2] 

    r1 = np.linspace(rp1,rc1, num = l, endpoint = 'true')
    r2 = np.linspace(rp2,rc2, num = l, endpoint = 'true')

    rho1 = eps*(1.0 + rp1**2)/(1.0 + r1**2) 

    rpd = np.zeros(len(time))
    for i in range(0, len(time)):
        rpd[i] = -0.5*(1.0/(np.sqrt(np.log(mt.cosh(1)))))*(mt.tanh(1 - time[i])/(np.sqrt(np.log(mt.cosh(1 - time[i])))))

    rpd1 = rpd[40]
    rpd2 = rpd[60]

    rdots, useless = rcdot(time,rc,rp)


    rcd1 = rdots[0,t1]

    u = ((2*rpd1*rp1)/(1.0 + rp1**2))*((1.0 + r1**2))*(np.arctan(r1) - 
       np.arctan(rc1)) + rcd1*(1.0 + r1**2)/(1.0 + rc1**2) 

    return u
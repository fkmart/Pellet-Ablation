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
    return rcd # out, r_pl 

"""
c = [-1.0, -5.0, -10.0, -100.0]
fig, ax = plt.subplots()
ax.set_xlabel(r'$\tilde{t}$', fontsize = 14)
ax.set_ylabel(r'$\frac{\mathrm{d}\tilde{A}}{\mathrm{d}\tilde{t}}$', fontsize = 14, rotation = 0)  
ax.yaxis.set_label_coords( -0.06, 0.45)
ax.xaxis.set_label_coords(0.50, -0.03)
plt.title('Modified ' + r'$D^2$' + ' law for varying ' + r'$c_2$')
for i in range(0,len(c)):
    dadt = dadt = c1*np.tanh((c[i]*(1-t)))
    plt.plot(t, dadt, label = r'$c_2 = ' + str(c[i])+'$')
    plt.legend()
plt.show()
"""
"""c2 = [-1,-5,-10,-100]
rplot = np.zeros((len(c2), len(t)))
rclot = np.zeros((len(c2), len(t)))
for i in range(0, len(c2)):
   print(i)
   rplot[i,:] = calc_rp(c2[i] ,t)
   rclot[i,:] = calc_rc(c2[i], t)"""
"""
rp = calc_rp(t)
rc = calc_rc(t)
rcd, rplot  = rcdot(t, rc, rp)

dt = t[1] - t[0]
fig, ax = plt.subplots()
ax.set_xlabel(r'$\tilde{r}$', fontsize = 14)
ax.set_ylabel(r'$\tilde{v}$', fontsize = 14, rotation = 0)
plt.title('Outflow velocity of cloud against cloud size for several times')

for i in range(10, len(t),20 ):
    plt.plot(rplot[i,:], rcd[i, :], label = r'$\tilde{t} = '+str(int(100*(dt*i))/100.0)+'$')
    
plt.legend()
plt.show()"""



"""fig, ax1  = plt.subplots()
ax1.set_xlabel(r'$\tilde{t}$', fontsize = 14)
ax1.set_ylabel(r'$\tilde{r}_p$', fontsize = 14, rotation = 0)
ax2 = ax1.twinx()
ax2.set_ylabel(r'$\tilde{r}_c$', fontsize = 14, rotation = 0)
linestyle = ['-', '--', '-.' , ':']
labels = []
for i in range(0,len(c2)):
    ax1.plot(t, rplot[i,:], 'royalblue', linestyle = linestyle[i])
    ax2.plot(t,rclot[i,:], 'forestgreen', linestyle = linestyle[i])
    plt.plot([],[], 'k', linestyle = linestyle[i], label=r'$c_2 =$' + str(c2[i]))
plt.legend(loc = 'center left')
plt.title('Pellet and Cloud radius as a function of time with changing parameter ' + str(r'$c_2$'))
plt.show()"""


"""
fig, ax = plt.subplots()
ax.set_xlabel(r'$\tilde{r}$', fontsize = 14)
ax.set_ylabel(r'$\tilde{\rho}$', fontsize = 14, rotation = 0)
ax.xaxis.set_label_coords(0.53, -0.04)
ax.yaxis.set_label_coords(-0.04, 0.54)
plt.title('Ablation cloud density against radius for several times')
tplot = np.linspace(0.0, 0.98, 50)
plt.yscale('log')
rp = calc_rp(tplot)
rc = calc_rc(tplot)
for i in range(5, len(tplot), 10):
    r = np.linspace(rp[i], rc[i], 500)
    den = eps*((1.0 + rp[i]**2)/ (1.0 + r[:]**2))
    ax.plot(r,den, label = r'$\tilde{t} = $' + str(round(tplot[i], 1)))
plt.legend()
plt.show()"""
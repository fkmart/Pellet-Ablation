import numpy as np 
import RK4 
import math as mt
import scipy.integrate as spit
import matplotlib.pyplot as plt

l = 1000
t = np.linspace(0, 1.0, l)
c1 = 1.0
c2 = -10.0

eps = 0.01
rc0 = 1.0

def rcloud(t,rc):
    rp = c1*(np.log(mt.cosh(c2*(1.0 - t)))/(np.log(mt.cosh(c2))))**(0.5)
    rpdot = -(c1*c2/(2.0))*((mt.tanh(c2*(1.0 - t))))*(np.log(mt.cosh(c2*(1 - t))))**(-0.5)
    rpdot *= (np.log(mt.cosh(c2)))**(-0.5)
    rcdot = ((eps -1)/eps)*((rp**2)/(rc**2))*rpdot*((1.0 + rc**2)/(1.0 + rp**2)) - 2*((rp*rpdot)/(rc**2))*((1.0 + rc**2)/(1.0 + rp**2))*(rc - rp + mt.atan(rp) - mt.atan(rc))
    return rcdot

rc = RK4.RK4_y_2(rcloud, rc0 , t)
rc0_arr = np.zeros(1) + 1.0
rc_auto = spit.solve_ivp(rcloud, (t[0],t[-2]), rc0_arr, method = 'RK45' )
print(rc[-1])
rc1 = rc_auto["y"]
rc1 = np.ravel(rc1)
print(rc1[-1])

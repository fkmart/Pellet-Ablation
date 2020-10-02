import numpy as np
from gen_var import t 
import matplotlib.pyplot as plt 
import scipy.integrate as spit

c1 = 1.0
c2 = 100.0
c3 = 1.0


def calc_rp(t):
    lt = len(t)
    dadt = c1*np.tanh((c2*(1-c3*t)))
    rnorm = c1*(np.log(np.cosh(c2*(1 - t)))/(np.log(np.cosh(c2))))**(0.5)
    drdtnorm = -0.5*c2*c1*((np.log(np.cosh(c2)))**(-0.5))*(((np.log(np.cosh(c2*(1 - t)))))**(-0.5))*(np.tanh(c2*(1-t)))
    np.save('rp_t'+str(lt)+'.npy', rnorm)
    return rnorm, dadt, drdtnorm 

rp, da, rpdot = calc_rp(t) 

k1 = 1.0
rc0 = np.asarray([1.0])
def calc_rc(time):
    lt = len(time)
    dt = time[1] - time[0]
    t = np.linspace(time[0], time[-1], len(time))
    dt = t[1] - t[0]
    eps = 0.01
    def func_rc(t,rc):
     
     drcdt = ((eps - 1.0)/eps)*((((c1*np.log(np.cosh(c2*(1 - t)))/(np.log(np.cosh(c2))))**(0.5))**2)/rc**2)*((1.0 
     # ((eps - 1.0)/eps)*((rnorm[i]**2)/(r**2))*((1.0 +
     + rc**2)/(1.0 +((c1*np.log(np.cosh(c2*(1 - t)))/(np.log(np.cosh(c2))))**(0.5))**2))*(-0.5*c2*c1)*((np.log(np.cosh(c2)))**(-0.5))*(((np.log(np.cosh(c2*(1.0 
     # r**2)/(1.0 + rnorm[i]**2))*drdtnorm[i]
     - t)))))**(-0.5))*(np.tanh(c2*(1-t))) - 2.0*((c1*(np.log(np.cosh(c2*(1.0 - t)))/(np.log(np.cosh(c2))))**(0.5))*-0.5*c2*c1*((np.log(np.cosh(c2)))**(-0.5))*(((np.log(np.cosh(c2*(1.0 
     - t)))))**(-0.5))*(np.tanh(c2*(1-t)))/rc**2)*((1.0 + rc**2)/(1.0 +c1*(np.log(np.cosh(c2*(1.0 - t)))/(np.log(np.cosh(c2))))**(1.0)))*(rc - 
     (np.log(np.cosh(c2*(1 - t)))/(np.log(np.cosh(c2))))**(0.5) + np.arctan((((np.log(np.cosh(c2*(1.0 
     - t)))/(np.log(np.cosh(c2))))**(0.5) -rc)/(1.0 +rc*(np.log(np.cosh(c2*(1 - t)))/(np.log(np.cosh(c2))))**(0.5)))))
     return drcdt 
    rcloud = spit.solve_ivp(func_rc, (t[0], t[len(t) -1]), rc0, method = 'RK45', max_step = dt)
    rc = rcloud["y"]
    #print(rc)
    rc = np.ravel(rc)
    rc = rc[:-1].copy()
    #rc = rc[::100]
    #print(rc)
    np.save('rc_t'+str(lt)+'.npy', rc)
    return rc
rc = calc_rc(t) 

def calc_rc_new(t, rp, rpdot):
    rc0 = np.asarray([1.0])
    dt = t[1] - t[0]
    eps = 0.01
    def rcfunc(t,rc):
        rcdot = (eps -1)*(rp*rp*rpdot/(rc**2))*((1.0 + rc**2)/(1.0 + rp**2)) - 2.0*((rp*rpdot/(rc*rc))*rc - rp + np.arctan((rp - rc/(1.0 +rp*rc))))
        return rcdot
    rcloud = spit.solve_ivp(rcfunc, (t[0], t[-1]), rc0, method = 'RK45', max_step = dt)
    rc = rcloud["y"]
    rc = np.ravel(rc)
    rc = rc[:-1].copy()
    return rc
#rc = calc_rc_new(t,rp,rpdot)
plt.plot(t, rp)
plt.plot(t, rc)
plt.show()

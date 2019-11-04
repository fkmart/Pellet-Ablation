import numpy as np 
import RK4 
import math as mt
import scipy.integrate as spit
import matplotlib.pyplot as plt
import os


l = 500
t = np.linspace(0, 1.0, l)
c1 = 1.0
#c2 = -10.0
c = [1.0, 5.0, 10.0, 100.0]
c = -np.asarray(c)
eps = 0.01
rc0 = 1.0 

title = 'rc'

path = os.path.join(os.path.expanduser('~'), 'Pictures')
filelist = [ f for f in os.listdir(path) if f.startswith(title)]
for f in filelist:
    os.remove(os.path.join(path, f))

def rcloud(t,rc):
    rp = c1*(np.log(mt.cosh(c2*(1.0 - t)))/(np.log(mt.cosh(c2))))**(0.5)
    rpdot = -(c1*c2/(2.0))*((mt.tanh(c2*(1.0 - t))))*(np.log(mt.cosh(c2*(1 - t))))**(-0.5)
    rpdot *= (np.log(mt.cosh(c2)))**(-0.5)
    rcdot = ((eps -1)/eps)*((rp**2)/(rc**2))*rpdot*((1.0 + rc**2)/(1.0 + rp**2)) - 2*((rp*rpdot)/(rc**2))*((1.0 + rc**2)/(1.0 + rp**2))*(rc - rp + mt.atan(rp) - mt.atan(rc))
    return rcdot

fig, ax = plt.subplots()
for i in range(0,len(c)):
    c2 = c[i]
    rc = RK4.RK4_y_2(rcloud, rc0 , t)
    ax.plot(t,rc, label = '$k_2 = $'  + str(-c[i]))
ax.set_xlabel(r'$\tilde{t}$', fontsize = 13)
ax.set_ylabel(r'$\tilde{r}_c$', fontsize = 13, rotation = 0.0)
ax.xaxis.set_label_coords(0.50, -0.04)
ax.yaxis.set_label_coords(-0.05, 0.46)
plt.legend()
fig.tight_layout()
plt.savefig(path+'/' + title + '.png', format = 'png', dpi = 1600)
plt.show()



#rc0_arr = np.zeros(1) + 1.0
#rc_auto = spit.solve_ivp(rcloud, (t[0],t[-2]), rc0_arr, method = 'RK45' )
print(rc[-1])
#rc1 = rc_auto["y"]
#rc1 = np.ravel(rc1)
#print(rc1[-1])

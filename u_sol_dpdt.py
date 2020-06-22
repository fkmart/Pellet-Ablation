from gen_var import rc,rp, eps
from gen_var import t as time 
import math as mt
import numpy as np
import matplotlib.pyplot as plt 
from stop_calc_rp_rc import rcdot

t1 = 20
t2 = 60 

shift = t2 - t1 

rp1,rp2 = rp[t1], rp[t2] 
rc1 ,rc2 = rc[t1], rc[t2]

l = 513
r1 = np.linspace(rp1,rc1, num = l, endpoint = 'true')
r2 = np.linspace(rp2,rc2, num = l, endpoint = 'true')

rho1 = eps*(1.0 + rp1**2)/(1.0 + r1**2)

rpd = np.zeros(len(time))
for i in range(0, len(time)):
    rpd[i] = -0.5*(1.0/(np.sqrt(np.log(mt.cosh(1)))))*(mt.tanh(1 - time[i])/(np.sqrt(np.log(mt.cosh(1 - time[i])))))

rpd1 = rpd[40]
rpd2 = rpd[60]

rdots, useless = rcdot(time,rc,rp)

rcd1 = rdots[t1,-1]

K1 = (1.0 + r1**2)*np.arctan(r1)
K2 = (1.0 + rc1**2)*np.arctan(rc1)
u = ((2*rpd1*rp1)/(1.0 + rp1**2))*((1.0 + r1**2)*np.arctan(r1) - (1.0 + rc1**2)*np.arctan(rc1)) + rcd1 

#revised u
u = ((2*rpd1*rp1)/(1.0 + rp1**2))*((1.0 + r1**2))*(np.arctan(r1) - np.arctan(rc1)) + rcd1*(1.0 + r1**2)/(1.0 + rc1**2)

fig,ax = plt.subplots()
ax.plot(r1, u)
plt.xlabel(r'$r$')
plt.ylabel(r'$u_r$', rotation = 0)
#plt.savefig('u_sol_t' + str(t1) + '.png', format = 'png', dpi = 1200)
plt.show()
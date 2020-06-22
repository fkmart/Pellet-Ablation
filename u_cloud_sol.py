import numpy as np 
import matplotlib.pyplot as plt 
import math as mt 
from stop_calc_rp_rc import rcdot 
from gen_var import rp, rc , eps, rpd, lt 
from gen_var import t as time

def u_sol(r,t1):
    import numpy
    from gen_var import rp,rc, rpd 
    rho = eps*(1.0 + rp[t1]**2)/(1.0 + r**2)

    rhodot = 2.0*rp[t1]*rpd[t1]/(1.0 + r**2) 

    rdots, useless = rcdot(time,rc,rp)
    rcd1 = rdots[t1,-1]

    u = rho[-1]*rcd1/rho + 2.0*(eps*rpd[t1]*rp[t1]/rho)*(np.arctan(r) - np.arctan(rc[t1]))
    return u 
"""
fig,ax = plt.subplots()
plt.xlabel(r'$\tilde{r}$')
plt.ylabel(r'$\tilde{u}_r$', rotation = 0)
for i in range(50,450,100): 
    r = np.linspace(rp[i], rc[i], num = 500, endpoint = 'point')
    u = u_sol(r,i)
    ax.plot(r, u, label = r'$\tilde{u}(\tilde{r},\tilde{t} = $'+'{:3.1f}'.format(i/lt) + ')')

plt.legend()
plt.savefig('u_sol.png', format = 'png', dpi = 1400)

plt.show()"""
import numpy as np 
import matplotlib.pyplot as plt 
from gen_var import rp, rc, eps, rpd, t 
from stop_calc_rp_rc import rcdot 
from u_cloud_sol import u_sol

t1 = 40 
rp1 = rp[t1] 
rc1 = rc[t1] 
rpd1 = rpd[t1] 

l = 500
r = np.linspace(rp1,rc1, num = l,endpoint = 'true')

alpha = 1e-7
rho = eps*(.0 + rp1**2)/(1.0 + r**2)
u_s = (1.0/rho)*(alpha - 2*eps*rpd1*rp1*(np.arctan(r) - np.arctan(rp)))
u_ns = u_sol(r,t1)

fig,ax = plt.subplots() 
ax.plot(r, u_ns)
ax.plot(r, u_s)
plt.show()
import numpy as np 
import scipy.special as scisp
import matplotlib.pyplot as plt 

l = 100
x = np.linspace(-1.0, 1.0, num = l, endpoint = 'true') # normalised to 1mm

n = np.zeros(l)
n[40:60] = 1.0 

#Transport with erfc/erf
D = 1e-3 #normalised to 1mm and 1ms
t=1 #normalised to 1ms

fig, ax = plt.subplots() 
dx = np.abs(x[60] - x[40])
#test case at t=0 
ts = 0 

ax.plot(x,n)
while (t<200):
    Z = x/(np.sqrt(4.0*D*t))
    Z_fac = dx/(np.sqrt(4.0*D*t))
    n_trans = 0.5*(scisp.erf(Z + Z_fac) - scisp.erf(Z - Z_fac))
    ax.plot(x,n_trans)
    t+=50
plt.show()
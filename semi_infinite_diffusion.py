import numpy as np
import matplotlib.pyplot as plt 
import scipy.special as scisp 
import romberg as ro

l = 257
x = np.linspace(0.0, 0.5, endpoint = 'true', num = l)

n = np.zeros(l)

n0 = 1.0 
n[0] =  n0 

""" In a semi-infinite plane, the diffusion solution is effectively just a exp(-x**2) solution.
There are additional constants detemrined by initial conditions such as amount of material, type
of material for diffusion constants etc """

B = (x[1] - x[0])*n0 *0.5
B = ro.romberg_samp(n,x)
D = 1e-5
t = 1

print('Initial check is ' + str(B))

fig, ax = plt.subplots()
ax.plot(x,n, label = r'$t = 0$' )
plt.xlim(0.0,0.1)
while t < 20:
    n_diff = (B/(np.sqrt(np.pi*D*t)))*np.exp(-(x**2)/(4.0*D*t))
    ax.plot(x,n_diff, label = r'$t = $' + str(t))
    check = ro.romberg_samp(n_diff,x)
    print('Check is ' + str(check) + 'at time ' + str(t))
    t += 5
plt.legend()
#plt.savefig('semi_inf_plane_diff.png', format = 'png', dpi = 1400)
plt.show()
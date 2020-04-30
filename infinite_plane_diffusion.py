import numpy as np 
import matplotlib.pyplot as plt 
import scipy.special as scisp 
import romberg as ro

l = 513
x = np.linspace(-1.0,1.0, num = l, endpoint = 'true')

n = np.zeros(l)
n0 = 1.0
n[int(l/2)] = n0

fig,ax = plt.subplots()
ax.plot(x,n, label = r'$t = 0$')

t = 5
D = 1e-5
B = (x[1] - x[0]) * n0
B = ro.romberg_samp(n,x)
print('Constant is ' + str(B))
while t < 200:
    Z = x/(np.sqrt(t*D*4.0))
    n_diff = 0.5*(B/(np.sqrt(np.pi*D*t)))*np.exp(-(x**2)/(4.0*D*t))
    ax.plot(x,n_diff, label = r'$t = $' + str(t))
    t += 50
    check = ro.romberg_samp(n_diff,x)
    print('Check is ' + str(check) + 'at time ' + str(t -50))
plt.legend()
plt.savefig('inf_plane_diff.png', format = 'png', dpi = 1400)
plt.show()
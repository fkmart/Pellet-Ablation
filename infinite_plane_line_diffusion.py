import numpy as np 
import matplotlib.pyplot as plt 
import scipy.special as scisp 
import scipy.signal as scisig
import romberg as ro

l = 257
x = np.linspace(-2.0,2.0, num = l, endpoint = 'true')
y = np.linspace(-1.0, 1.0, num = l, endpoint = 'true')
n = np.zeros(l)

l1 = 33
n[110:110 + l1] = np.linspace(0.5,2.5, num = l1)
#n_new = np.linspace(0.5,1.2, num = l1)

#y = np.copy(x[110:110 + l1])
fig,ax = plt.subplots()
ax.plot(x,n, label = r'$t = 0$')

t = 1
D = 1e-5

B = ro.romberg_samp(n,x)
n_diff = np.zeros(l)
print('Constant is ' + str(B))
while t < 200:
    n *= 1.0
    Z = x/(np.sqrt(t*D*4.0))
    f = (np.exp(-Z**2))/(np.sqrt(np.pi*D*t))
    n_diff = scisig.convolve(n, f, mode = 'same')
    j = np.sum(f)
    print('sum equals ' + str(j))
    ax.plot(x,B*n_diff/j, label = r'$t = $' + str(t))
    t += 50
    check = ro.romberg_samp(n_diff/j,x)
    print('Check is ' + str(check) + ' at time ' + str(t -50))
plt.legend()
plt.show()
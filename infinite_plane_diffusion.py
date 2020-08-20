import numpy as np 
import matplotlib.pyplot as plt 
import scipy.special as scisp 
import romberg as ro
import scipy.signal as scisig

l = 513
x = np.linspace(-0.5,0.5, num = l, endpoint = 'true')
y = np.linspace(-0.1,0.1, num = l, endpoint = 'true')
n_new = np.ones(len(y))
n = np.zeros(l)
n0 = 1.0
xl = int(l/2) - 20 
xr = xl + 40
n[xl:xr] = n0

fig,ax = plt.subplots()
ax.plot(x,n, label = r'$\tilde{n}(t = 0)$')
ax.set_xlabel(r'$\tilde{x}$', rotation = 0, fontsize = 12)
ax.set_ylabel(r'$\tilde{n}$', rotation = 0, fontsize = 12)
ax.xaxis. set_label_coords(0.45,-0.04)
ax.yaxis.set_label_coords(-0.04, 0.48)
t = 50
D = 1e-5
B = ro.romberg_samp(n,x)
print('Constant is ' + str(B))
c = 1
n_diff_out = np.copy(n)
dx = x[1] - x[0]
while c < 5:
    Z = x/(np.sqrt(t*D*4.0))
    n_diff = (1.0/(np.sqrt(4.0*np.pi*D*t)))*np.exp(-Z**2)
    #n_diff_out = scisig.convolve(n_diff,n, mode = 'same', method = 'direct')
    n_diff_out = dx*np.convolve(n_diff,n_diff_out, mode = 'same')
    check = ro.romberg_samp(n_diff_out,x)
        
    #n_diff_out *= B/check
    final_check = ro.romberg_samp(n_diff_out,x)
    ax.plot(x,n_diff_out, label = r'$\tilde{n}(t = $'+str(c) + r'$\Delta\tilde{t})$')
    #t += 50
    c +=1
    print('Check is ' + str(final_check) + ' at time ' + str(t -50))
plt.legend()
#plt.savefig('inf_plane_diff.png', format = 'png', dpi = 1400)
plt.show()
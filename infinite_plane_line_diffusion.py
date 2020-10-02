import numpy as np 
import matplotlib.pyplot as plt 
import scipy.signal as scisig
import romberg as ro
from gen_var import m_h, m_h2, A_con, omega, s 

l = 513
x = np.linspace(-3.0, 8.0, num = l, endpoint = 'true')

n = np.zeros(l)

l1 = 60
il = 270
iu = il+l1
n[il:iu] = np.linspace(0.5,2.5, num = l1)
#n_new = np.linspace(0.5,1.2, num = l1)

#y = np.copy(x[110:110 + l1])
fig,ax = plt.subplots()
plt.xlabel(r'$\tilde{x}$', fontsize = 12)
ax.xaxis.set_label_coords(0.45,-0.04)
plt.ylabel(r'$\tilde{n}$', fontsize = 12, rotation = 0)
ax.yaxis.set_label_coords(-0.06,0.45)
ax.plot(x,n, label = 'initial')

t = 5*10**(-4)

T = 11600
#D = 1e-5

D = A_con*T**(1.5) * np.sqrt(m_h2**-1 + m_h**-1)/(1.0*s**2 *omega)
D *= 10**(-4)
print(D)
B = ro.romberg_samp(n,x) # mass in system
n_diff = np.zeros(l)
print('Constant is ' + str(B))
c = 1
r_cent = (x[iu] + x[il])*0.5
r_cent = x[0]+ 0.5*(-x[0] + x[-1])
while t < 5*10**(-3):
    
    Z = (x-r_cent)/(np.sqrt(t*D*4.0))
    f = (np.exp(-Z**2))/(np.sqrt(4.0*np.pi*D*t))
    n_diff = scisig.convolve(n, f, mode = 'same')
    n_diff *= (x[1] - x[0])
    j = ro.romberg_samp(n_diff,x)
    print('Mass after diffusion = ' + str(j))
    ax.plot(x,n_diff, label = 't = ' + str(c))
    t += 10**(-3)
    check = ro.romberg_samp(n_diff,x)
    print('Check is ' + str(check) + ' at time ' + str(t -5*10**(-5)))
    c +=1
plt.legend()
#plt.savefig('diff_over_t.png', format = 'png', dpi = 1400)
plt.show()
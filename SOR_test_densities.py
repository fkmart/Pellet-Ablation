import numpy as np 
import iterative_sol as SOR
import discret_mat
from gen_var import rc, rp  , r0, e , n_r, r, phi_p, epsilon0
import matplotlib.pyplot as plt


savedir_an = '/home/kyle/Documents/Python_Code/stop_code/one_size_fits_all/one_iteration/analysed_outputs/'

i = 50
dens = -np.loadtxt(savedir_an + 'real_density_pot_test_t'+str(i) +'pot0.0'+'.txt')

neut_dens = np.zeros(n_r)
lr = len(r)
low = next(p[0] for p in enumerate(r) if p[1] > rp[i])
up = next(p[0] for p in enumerate(r) if p[1] > rc[i])

l = len(dens[:,0])

r_dom = r[up - l : up]
#r_dom[:] = r_dom[:] - r_dom[0]
#r_dom[:] /= r_dom[-1]
A = discret_mat.discret(r_dom) #r_domain is from just in front of the pellet to the cloud at the new time.
norm = r0**2*(e*10**19)/(10**3 * epsilon0)
pot_in = np.zeros(l)
pot_in[-1] = 0.0
pot_in[0] = 0.0
phic = SOR.SOR(A, pot_in, dens[:,0], r_dom)

fig, ax  = plt.subplots() 
ax.plot(r_dom, phic*10**3, color = 'royalblue', label = r'$\phi_c$')
ax.set_ylabel(r'$\phi_c$', fontsize = 12, rotation = 0.0)
ax.set_xlabel(r'$\tilde{r}$', fontsize = 12)
ax2 = ax.twinx() 
ax2.plot(r_dom, dens[:,0], color = 'chocolate', label = r'$\tilde{n}$')
ax2.set_ylabel(r'$\tilde{n}_c$', fontsize = 12, rotation = 0.0)
plt.show()
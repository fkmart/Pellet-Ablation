import numpy as np 
import matplotlib.pyplot as plt 
import os 
from gen_var import rp, rc, pel_pot, lp, sig, t_start
import gauss_test_pot as gtp

t_start = 50
direc = os.getcwd()
#sub_dir = '/one_iteration_phic/analysed_outputs/'
sub_dir = os.path.join(direc, 'one_iteration_phic', 'analysed_outputs') + os.sep
r_int = np.linspace(rp[t_start], rc[t_start], 602)
savedir = os.path.join(direc,'pictures') + os.sep

mid = 301
ind1 = 15
ind = 150

fig, ax  = plt.subplots(figsize = (10,8))

c = ['forestgreen','midnightblue', 'saddlebrown', 'y', 'darkred', 'dodgerblue', 'peru', 'gold', 'limegreen', 'tomato', 'firebrick']
y = str(1.00)
plt.axvline(r_int[mid], color = 'black', linestyle = '--')
file = np.loadtxt(sub_dir + 'density_pot_test_t50pot-0.0sig1.0.txt')
ax.plot(file[:,0], file[:,1], color = 'black', label = 'Pure CSDA')
u = 0
for p in range(2, 12,2):
    pot = gtp.gauss_func(pel_pot[p], float(y), r_int[mid], r_int)
    x = format(pel_pot[p], '.1f') 
    file = np.loadtxt(sub_dir + 'density_pot_test_t' + str(t_start) + 'pot' + x + 'sig' + y + '.txt')
    ax.plot(file[:,0], file[:,1], color = c[u], label = r'$\phi_{\mathrm{max}} = $' + str(int(pel_pot[p])) + 'V')
    u +=1
plt.legend(ncol = 2, loc = 'lower right')

ax.set_xlabel(r'$\tilde{r}$', fontsize = 12)
ax.xaxis.set_label_coords(0.51,-0.03)
ax.set_ylabel(r'$\tilde{n}$', rotation = 0, fontsize = 12)
#ax.set_ylabel(r'$n_e / \ 10^{19}\mathrm{m}^{-3}$',fontsize = 10, rotation = 0)
ax.yaxis.set_label_coords(-0.06, 0.47)
ax.set_yscale('log')

ax.set_yscale('log')
#ax2 = ax.twinx()
#ax2.set_ylabel(r'$\phi$/V',fontsize = 12, rotation = 0)
#ax2.yaxis.set_label_coords(1.03, 0.55)
#for p in range(0, lp):
    #pot = gtp.gauss_func(pel_pot[p],1.00, r_int[mid], r_int)
    #ax2.plot(r_int[-ind:], pot[-ind:], color = c[p], linestyle = '--')

plt.text(0.8,0.0006, r'$\leftarrow$' + 'to pellet')
plt.text(5.5, 0.0006, r'$\rightarrow$' + 'to plasma')

y = float(y)
y = str(int(y))
plt.savefig(savedir + 'bins_pot_sig.png', format = 'png', dpi = 1200)
#plt.grid('on')
plt.show()
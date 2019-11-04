import numpy as np 
import matplotlib.pyplot as plt
import os
from gen_var import rp,rc,t,eps, dt

fig, ax = plt.subplots() 
div = int(0.2*len(t))
i_start = int(0.3*len(t))
i_end = len(t)

l = 500

title = 'neut_dens'

mydir = '~/Pictures'
path = os.path.join(os.path.expanduser('~'), 'Pictures')
filelist = [ f for f in os.listdir(path) if f.startswith(title)]
for f in filelist:
    os.remove(os.path.join(path, f))

for i in range(i_start, i_end, div):
    r = np.linspace(rp[i], rc[i], l)
    dens = eps*((1.0 + rp[i]**2)/(1.0 + r**2))
    ax.plot(r,dens,label = r'$\tilde{t} = $'+str(round(dt*i,1)))
ax.set_xlabel(r'$\tilde{r}$', fontsize = 13)
ax.set_ylabel(r'$\tilde{\rho}_c$', fontsize = 13, rotation = 0)
ax.xaxis.set_label_coords(0.50, -0.04)
ax.yaxis.set_label_coords(-0.04, 0.50)
ax.set_yscale('log')
plt.legend()
fig.tight_layout()
plt.show()
plt.savefig(path+'/' + title + '.eps', format = 'eps', dpi = 1600)



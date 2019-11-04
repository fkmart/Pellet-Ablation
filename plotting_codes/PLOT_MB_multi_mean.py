import numpy as np  
import matplotlib.pyplot as plt 
import os 

title = 'multi_MB'

save_path = os.path.join(os.path.expanduser('~'), 'Documents/thesis/Local_figures/')
filelist = [ f for f in os.listdir(save_path) if f.startswith(title)]
for f in filelist:
    os.remove(os.path.join(save_path, f))


E = np.arange(10,10000,10)
eb = [500,750, 1000,3000]
fig, ax  = plt.subplots()

for i in range(0, len(eb)):
    e_bar = eb[i]
    MB_ener = ((2.07*E**(0.5))/e_bar**(1.5))*np.exp(-1.5*(E/e_bar))
    ax.plot(E, MB_ener, label = r'$\bar{E} = $'+ str(e_bar) + 'eV')
#ax.set_yscale('log')
ax.set_xlabel(r'$E/\mathrm{eV}$', fontsize = 12)
ax.set_ylabel(r'$f_{MB}$', fontsize = 12, rotation = 0)
ax.xaxis.set_label_coords(0.5, -0.04)
ax.yaxis.set_label_coords(-0.06, 0.47)
plt.legend()
plt.savefig(save_path + title + '.png', format = 'png', dpi = 1200)
plt.show()
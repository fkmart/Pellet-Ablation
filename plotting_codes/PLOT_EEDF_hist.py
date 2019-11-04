import numpy as np
import matplotlib.pyplot as plt 
import electron 
import os

savedir = os.path.join(os.path.expanduser('~'), 'Documents/thesis/Local_figures')
res = 50
ener = np.arange(0,7000,res)
mid, bins, MB_norm = electron.dist_calc(ener,res,1000)

fig, ax = plt.subplots()
ax2 = ax.twinx()
ax.bar(mid,bins, width = res, edgecolor = 'black')
ax2.set_ylim(bottom = 0.0, top = 0.00075)
ax2.plot(ener,MB_norm, color = 'firebrick')
ax.set_ylabel(r'$\int_E ^ {E + \Delta E} f_{MB}(E) \mathrm{d}E$', fontsize = 7, rotation = 0)
ax.yaxis.set_label_coords(-0.08, 0.44)
ax2.set_ylabel(r'$f_{MB}(E)$', rotation = 0, fontsize = 10)
ax2.yaxis.set_label_coords(1.06, 0.5)
ax.set_xlabel(r'$E/\mathrm{eV}$', fontsize = 12)
plt.savefig(savedir+'/eedf_hist.png', format = 'png', dpi = 1200)

plt.show()
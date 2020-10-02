# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 15:44:23 2020

@author: Kyle
"""
import numpy as np
import matplotlib.pyplot as plt 
import os 
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 31 11:16:35 2020

@author: Kyle
"""
import os 
import matplotlib.pyplot as plt 
import numpy as np 
from gen_var import delta_t, RME, M_fac

direc = os.getcwd() 
TYPE = 'ion'

l = 2000
step = 400
fig, ax  = plt.subplots()
ax.set_xlabel(r'$\tilde{r}$', fontsize = 12, rotation = 0 )
ax.set_ylabel(r'$\left|\phi\right|$/V', fontsize = 12, rotation = 0)
ax.xaxis.set_label_coords(0.45,-0.03)
ax.yaxis.set_label_coords(-0.05,0.52)

for j in range(1, l, step):
    load_dir = os.path.join(direc, 'many_iteration', TYPE, '1000eV', 'analysed_outputs') + os.sep
    file = np.genfromtxt(load_dir + 'outputs' + str(j) + '.txt', delimiter = ',', dtype = 'str')
    pot = file[5:,3]
    pot = np.asarray([float(w) for w in pot])
    r = file[5:,0]
    r = np.asarray([float(w) for w in r])
    ax.plot(r,np.abs(pot)*RME*M_fac, label = r'$\Delta \tilde{t} = $' +'{:3.2f}'.format((j-1)*delta_t))
plt.yscale('log')
plt.xlim(r[0], 6.5)
#plt.xscale('log')
plt.legend(ncol = 2, loc = 'lower center')
plt.savefig('potential_' + TYPE +'_log.png', format = 'png', dpi = 1400)
plt.show()
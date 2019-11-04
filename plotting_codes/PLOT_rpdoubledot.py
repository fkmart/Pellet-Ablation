import numpy as np 
import matplotlib.pyplot as plt 
import math as mt 
import os 

title = 'rpddot'
path = os.path.join(os.path.expanduser('~'), 'Documents/thesis/Local_figures/')
filelist = [ f for f in os.listdir(path) if f.startswith(title)]
for f in filelist:
    os.remove(os.path.join(path, f))


l = 5000
k = [1.0, 5.0, 10.0, 100.0]
k = np.asarray(k)
t = np.linspace(0.0, 1.0, l) 

x = np.zeros(len(t))
rpddot = np.zeros(len(t))
fig, ax = plt.subplots()
for i in range(0, len(k)):
    for j in range(0, len(t)):
        x[j] = k[i]*(1.0 - t[j])
        rpddot[j] = -0.5*k[i]*(1.0/(np.sqrt(np.log(mt.cosh(k[i])))))
        rpddot[j] *= ((np.log(mt.cosh(x[j])))**(-0.5))*(1.0/(mt.cosh(x[j])*mt.cosh(x[j])))*(-k[i]) + mt.tanh(x[j])*(0.5*k[i])*((np.log(mt.cosh(x[j])))**(-1.5))*mt.tanh(x[j])
    ax.plot(t, rpddot, label = r'$k_2 = $'+ str(k[i]))
ax.set_xlabel(r'$\tilde{t}$', fontsize = 12)
ax.set_ylabel(r'$\tilde{\ddot{r}}_p$', rotation = 0, fontsize = 12)
ax.xaxis.set_label_coords(0.5, -0.03)
plt.legend()
plt.savefig(path + title + '.png', format = 'png', dpi = 1200)
plt.show()
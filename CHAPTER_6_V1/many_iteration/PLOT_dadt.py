import numpy as np 
import matplotlib.pyplot as plt 
import math as mt
import os

l = 500
t = np.linspace(0.0,1.0, l)
k = [1.0,5.0,10.0,100.0]
k = -np.asarray(k)

title = 'dadt'

mydpi = 1600

path = os.path.join(os.path.expanduser('~'), 'Pictures')
filelist = [ f for f in os.listdir(path) if f.startswith(title)]
for f in filelist:
    os.remove(os.path.join(path, f))

fig, ax = plt.subplots()
#fig(figsize = (1.600,2.400), dpi = 100)
for i in range(0,len(k)):
    dadt = np.tanh(k[i]*(1.0 - t))
    ax.plot(t,dadt, label = r'$k_2 = $' + str(k[i]))
ax.set_xlabel(r'$\tilde{t}$', fontsize = 13)
ax.set_ylabel(r'$\frac{\mathrm{d}\tilde{A}}{\mathrm{d}\tilde{t}}$', fontsize = 14, rotation = 0.0)
ax.xaxis.set_label_coords(0.5,-0.04)
ax.yaxis.set_label_coords(-0.06, 0.46)
plt.legend()
fig.tight_layout()

plt.show()
#plt.savefig(path + '/' + title + '.png', format = 'png', dpi = mydpi)
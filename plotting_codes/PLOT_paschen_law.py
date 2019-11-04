import numpy as np
import matplotlib.pyplot as plt


C1 = 7.95
C2 = 263

#NEED TO LEARN WHAT GAMMA IS FOR H2
gamma= 1
pd = np.linspace(0.1, 10, 500)


Vb = C2*pd/(np.log((C1*pd) - np.log(1.0 + 1.0/gamma)))

fig, ax = plt.subplots()

ax.plot(pd,Vb)
ax.set_xscale('log')
plt.show() 
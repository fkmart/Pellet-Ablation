import numpy as np 
import matplotlib.pyplot as plt 

x = np.linspace(0.0, 10.0, 10)
y = x**2 

fig, ax = plt.subplots() 
ax.plot(x,y)
plt.savefig('dummy.png', format = 'png', dpi = 800)
plt.show()
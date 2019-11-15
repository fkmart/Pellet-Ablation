import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits import mplot3d

fig = plt.figure()
ax = plt.axes(projection = '3d') 
x = np.linspace(0, 4*np.pi,100)
z = np.linspace(0,5,5)
y = np.zeros((len(x),len(z))) 
for i in range(0, len(z)):
    y[:,i] = np.sin(x[:]*z[i])
    ax.plot3D(x,y[:,i], z[i])
ax.view_init(20,170)
plt.show()

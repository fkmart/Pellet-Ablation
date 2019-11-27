import numpy as np 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D 

x = np.linspace(-5.0, 5.0, endpoint = 'true', num = 50)
y = np.linspace(-5.0, 5.0, endpoint = 'true', num = 50)

X,Y = np.meshgrid(x,y)
Z = np.sqrt(X**2 + Y**2 )
Z = np.sin(Z)
C = np.linspace(0.0,10.0, num = Z.size)
C = np.reshape(C, np.shape(Z), 'F')
scamap = plt.cm.ScalarMappable(cmap = 'viridis')
colors = scamap.to_rgba(C)
fig = plt.figure(figsize = (10.0,8.0))
ax= plt.axes(projection = '3d')
ax.plot_surface(X,Y,Z, facecolors = colors, cmap = 'viridis')
fig.colorbar(scamap)
plt.show()

import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D 
import matplotlib.tri as mtri

#ellipse params 

a = 1.0
b = 2.0

x = np.linspace(-a, a, endpoint = 'true', num = 51)
y = np.linspace(-b, b, endpoint = 'true', num = 51)

X,Y = np.meshgrid(x,y)
#X,Y = X.flatten(), Y.flatten()
z = (X/a)**2 + (Y/b)**2

"""
#Define the colourmesh 
tri = mtri.Triangulation(X, Y)
fig = plt.figure()
plt.figaspect(a/b)
ax = fig.add_subplot(1, 2, 1, projection='3d')
ax.plot_trisurf(X, Y, z, triangles=tri.triangles, cmap=plt.cm.Spectral)
el = b*np.sqrt(1.0 - (x/a)**2)
ax.plot(x,el,np.ones(len(x)), color = 'black', linestyle = '--')
plt.show()"""

fig, ax  = plt.subplots() 
ax = ax.pcolormesh(X,Y,z, cmap = 'viridis')
plt.show()

import numpy as np 

def rk(x,y,func,h):
    k1 = np.abs(h)*func(x,y) #next four lines are RK4 
    k2 = np.abs(h)*func(x + h*0.5, y + k1*0.5)
    k3 = np.abs(h)*func(x + h*0.5, y + k2*0.5)
    k4 = np.abs(h)*func(x + h, y + k3)
    y_new = y + (1.0/6.0)*(k1 + 2.0*k2 + 2.0*k3 + k4)
    return y_new, k1,k2,k3,k4
import numpy as np

def RK4_t(f,x0,h,y_end, y0):
    y = y0
    y_out = []
    y_out.append(y0)
    x = [x0]
    while (y > y_end):
        k1 = h*f(y)
        k2 = h*f(y+0.5*k1)
        k3 = h*f(y+0.5*k2)
        k4 = h*f(y+k3)
        y = (1.0/6.0)*(k1+2*k2 + 2*k3 + k4)
        y_out.append(y)
        x.append(x[-1]+h)
    return y_out, x

def RK4_y(f, y0, t):
    y_out = []
    y_out.append(y0)
    c = 1
    h = t[1] - t[0]
    y = y0
    
    while (c < len(t)):
        h = t[c] - t[c-1]
        i = t[c]
        k1 = h*f(i)
        k2 = h*f(i+0.5*k1)
        k3 = h*f(i+0.5*k2)
        k4 = h*f(i + k3)
        y = y + (1.0/6.0)* (k1 + 2*k2 + 2*k3 + k4)
        y_out.append(y)
        c+=1
    y_out = np.asarray(y_out)
    return y_out

def RK4_y_2(f, y0, t):
    y_out = []
    y_out.append(y0)
    c = 1
    h = t[1] - t[0]
    y = y0
    #t = t[:-2]
    
    while (c < len(t)):
        h = t[c] - t[c-1]
        i = t[c]
        k1 = h*f(i,y)
        k2 = h*f(i+ 0.5*h,y+0.5*k1)
        k3 = h*f(i + 0.5*h,y+0.5*k2)
        k4 = h*f(i + h, y + k3)
        y = y + (1.0/6.0)* (k1 + 2*k2 + 2*k3 + k4)
        y_out.append(y)
        
        c+=1
    y_out = np.asarray(y_out)
    return y_out
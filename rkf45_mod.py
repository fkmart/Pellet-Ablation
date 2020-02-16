import numpy as np

def rkf(x,y,h,func,xr,xl):
    # min and max time step
    hmin = 1e-5
    hmax = 5e-1
    # min and max errors
    err_min = 1e-7
    err_max = 1e-5
    # max number of iterations
    N = int((xr - xl)/hmin)
    N = 100
    if x+h > xr:
        h = xr-x
    #butcher-tableau 
    BT = np.zeros((6,6))
    BT[1,:] = [0.25,0.25,0.0, 0.0, 0.0, 0.0]
    BT[2,:] = [0.125, 3.0/32.0, 9.0/32.0, 0.0, 0.0, 0.0]
    BT[3,:] = [12.0/13.0, 1932.0/2197.0, -7200.0/2197.0, 7296.0/2197.0, 0.0, 0.0]
    BT[4,:] = [1.0, 439.0/216.0, -8.0, 3680.0/513.0, -845.0/4104.0, 0.0]
    BT[5,:] = [0.5, -8.0/27.0, 2.0, -3544.0/2565.0, 1859./4140.0 , -11.0/40.0]
    
    #butcher-tableau for solutions

    BTS = np.zeros((2,6))
    BTS[0,:] = [25.0/216.0, 0.0, 1408.0/2565.0, 2197.0/4104.0, - 0.2, 0.0]
    BTS[1,:] = [16.0/135.0, 0.0, 6656.0/12825.0,28561.0/56430.0, -9.0/50.0, 2.0/55.0 ]

    #calculate the k's 
    for i in range(N):
        k1 = func(x,y)
        k2 = func(x + BT[1,0]*h, y + h*(k1*BT[1,1]))
        k3 = func(x + BT[2,0]*h, y + h*(k1*BT[2,1] + k2*BT[2,2]))
        k4 = func(x + BT[3,0]*h, y + h*(k1*BT[3,1] + k2*BT[3,2] + k3*BT[3,3]))
        k5 = func(x + BT[4,0]*h, y + h*(k1*BT[4,1] + k2*BT[4,2] + k3*BT[4,3] + k4*BT[4,4]))
        k6 = func(x + BT[5,0]*h, y + h*(k1*BT[5,1] + k2*BT[5,2] + k3*BT[5,3] + k4*BT[5,4] + k5*BT[5,5]))

        #calculate fourth and fifth order solutions 
        y4 = y + h*(BTS[0,0]*k1 + BTS[0,2]*k3 + BTS[0,3]*k4 + BTS[0,4]*k5)
        y5 = y + h*(BTS[1,0]*k1 + BTS[1,2]*k3 + BTS[1,3]*k4 + BTS[1,4]*k5 + BTS[1,5]*k6)

        err = np.abs(y4-y5)

        if err < err_min:
                # if error small, enlarge h, but match final simulation time
            h = min(2.*h,hmax)            
            if x+h > xr:
                h = xr-x
                break
        elif err > err_max:
                # if error big, reduce h
                h = max(h/2.,hmin)
        else:
                # error is ok, take this h and y5
                break
        
    if i==N-1:
        print("max number of iterations reached, check parameters")
            
    return x+h, y5, h , err

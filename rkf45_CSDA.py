import numpy as np
from electron import RME, M_fac 
from gen_var import I, trunc_fac
import frigg_it_logic as fil
import step_calc as sc
import sys
def rkf(x,y,h,func,BT,BTS,xr,xl,I_non):
    # min and max time step
    hmin = 1e-6
    hmax = 5e-2
    h_arr = []
    br_no = 20
    # min and max errors
    if y > 10000.0/(RME*M_fac): #1e-8 and 1e-4 works
        err_min = 1e-10
        err_max = 1e-4

    elif y > 1000.0/(RME*M_fac): # 1e-10 and 1e-6 works - not normalised
        err_min = 1e-10
        err_max = 1e-6

    else: # 1e-13 and 1e-8 works - not normalised
        err_min = 1e-12
        err_max = 1e-7

    #err_min = 1e-11
    #err_max = 1e-4
    # max number of iterations
    N = int((xr - xl)/hmin)
    N = 100
    #if x+h > xr:
    #    h = xr-x
    if x - np.abs(h) < xl:
        h = x - xl 
    #calculate the k's 
    h = -np.abs(h)
    k1 = func(x,y)
    k2 = func(x + BT[1,0]*h, y - h*(k1*BT[1,1])) # should be + not - but need to account for negative dr 
    k3 = func(x + BT[2,0]*h, y - h*(k1*BT[2,1] + k2*BT[2,2]))
    k4 = func(x + BT[3,0]*h, y - h*(k1*BT[3,1] + k2*BT[3,2] + k3*BT[3,3]))
    k5 = func(x + BT[4,0]*h, y - h*(k1*BT[4,1] + k2*BT[4,2] + k3*BT[4,3] + k4*BT[4,4]))
    k6 = func(x + BT[5,0]*h, y - h*(k1*BT[5,1] + k2*BT[5,2] + k3*BT[5,3] + k4*BT[5,4] + k5*BT[5,5]))
    k = [k1,k2,k3,k4,k5,k6]
    for i in range(N):

        #calculate fourth and fifth order solutions 
        y4 = y - h*(BTS[0,0]*k[0] + BTS[0,2]*k[2] + BTS[0,3]*k[3] + BTS[0,4]*k[4])
        y5 = y - h*(BTS[1,0]*k[0] + BTS[1,2]*k[2] + BTS[1,3]*k[3] + BTS[1,4]*k[4] + BTS[1,5]*k[5])

        #y_norm = np.max(np.asarray(y4,y5))
        #err = np.abs(y4/y_norm-y5/y_norm) # need to remove y norms if this doesn't work
        err = np.abs(y4 - y5)

        if err < err_min and y4 > 0.0 and y5 > 0.0:
                # if error small, enlarge h, but match final simulation time
            h = min(2.*np.abs(h),hmax)            
            if x - h < xl:
                h = xl - x
                break
        elif err > err_max or  y4 < 0.0 or y5 < 0.0:
                # if error big or energy is less than critical value, reduce h
                h = max(np.abs(h/2.),hmin)
        else:
                # error is ok, take this h and y5
                break
        h_arr.append(h)
        if i%br_no ==0 and i>0:
            test, h_test = fil.frigg_it_logic(h_arr, br_no)
            if test == 'true':
                h = h_test
                y4 = y - h*(BTS[0,0]*k1 + BTS[0,2]*k3 + BTS[0,3]*k4 + BTS[0,4]*k5)
                y5 = y - h*(BTS[1,0]*k1 + BTS[1,2]*k3 + BTS[1,3]*k4 + BTS[1,4]*k5 + BTS[1,5]*k6)
                err = np.abs(y4-y5)
                #print("'Frigg it' logic used")
                break
            else:
                pass
    
        h = -np.abs(h)
        y_steps = np.ones(5)
        subvert = 'false'
        for p in range(0,5):
            y_steps[p] = sc.ys_calc(h,k,BT,p)
            if y - y_steps[p] < trunc_fac*I_non:
                h *= 0.5
                subvert = 'true'
                break
            else:
                k[p] = sc.k_calcs(x,y,h,BT,func,p,k)
        if subvert =='true':
            break
        else:
            pass
    for p in range(0, 6):
        k[p] = sc.k_calcs(x,y,h,BT,func,p,k)
    y4 = y - h*(BTS[0,0]*k[0] + BTS[0,2]*k[2] + BTS[0,3]*k[3] + BTS[0,4]*k[4])
    y5 = y - h*(BTS[1,0]*k[0] + BTS[1,2]*k[2] + BTS[1,3]*k[3] + BTS[1,4]*k[4] + BTS[1,5]*k[5])
    err = np.abs(y4-y5)
    if i==N-1:
        #print("max number of iterations reached, check parameters, using h_min")
        h = -hmin
        k1 = func(x,y)
        k2 = func(x + BT[1,0]*h, y - h*(k1*BT[1,1])) # should be + not - but need to account for negative dr 
        k3 = func(x + BT[2,0]*h, y - h*(k1*BT[2,1] + k2*BT[2,2]))
        k4 = func(x + BT[3,0]*h, y - h*(k1*BT[3,1] + k2*BT[3,2] + k3*BT[3,3]))
        k5 = func(x + BT[4,0]*h, y - h*(k1*BT[4,1] + k2*BT[4,2] + k3*BT[4,3] + k4*BT[4,4]))
        k6 = func(x + BT[5,0]*h, y - h*(k1*BT[5,1] + k2*BT[5,2] + k3*BT[5,3] + k4*BT[5,4] + k5*BT[5,5]))
        k = [k1,k2,k3,k4,k5,k6]
        y4 = y - h*(BTS[0,0]*k[0] + BTS[0,2]*k[2] + BTS[0,3]*k[3] + BTS[0,4]*k[4])
        y5 = y - h*(BTS[1,0]*k[0] + BTS[1,2]*k[2] + BTS[1,3]*k[3] + BTS[1,4]*k[4] + BTS[1,5]*k[5])

        #y_norm = np.max(np.asarray(y4,y5))
        #err = np.abs(y4/y_norm-y5/y_norm) # need to remove y norms if this doesn't work
        err = np.abs(y4 - y5)
        #sys.exit()
    #k = [k1,k2,k3,k4,k5,k6]
    err = np.abs(y4 - y5)
    return x - np.abs(h), y5, h , err, k
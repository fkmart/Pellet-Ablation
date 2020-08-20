import numpy as np 
from gen_var import r0cgs, solid_dens, zovera,I, eps
from electron import RME, M_fac
from numba import jit, njit

def k_calcs(x,y,h,BT, j, k, func) :
        
    if j ==0 :
        k_out = func(x,y)
    elif j==1:
        k_out = func(x + BT[1,0]*h, y - h*(k[0]*BT[1,1])) # should be + not - but need to account for negative dr 
    elif j==2:
        k_out = func(x + BT[2,0]*h, y - h*(k[0]*BT[2,1] + k[1]*BT[2,2]))
    elif j==3:
        k_out = func(x + BT[3,0]*h, y - h*(k[0]*BT[3,1] + k[1]*BT[3,2] + k[2]*BT[3,3]))
    elif j==4:
        k_out = func(x + BT[4,0]*h, y - h*(k[0]*BT[4,1] + k[1]*BT[4,2] + k[2]*BT[4,3] + k[3]*BT[4,4]))
    elif j==5:
        k_out = func(x + BT[5,0]*h, y - h*(k[0]*BT[5,1] + k[1]*BT[5,2] + k[2]*BT[5,3] + k[3]*BT[5,4] + k[4]*BT[5,5])) 
    return k_out 

def ys_calc(h,k,BT,j):
    if j == 0:
        ys = h*(k[0]*BT[1,1])
    elif j ==1:
        ys = h*(k[0]*BT[2,1] + k[1]*BT[2,2])
    elif j==2:
        ys = h*(k[0]*BT[3,1] + k[1]*BT[3,2] + k[2]*BT[3,3])
    elif j ==3: 
        ys = h*(k[0]*BT[4,1] + k[1]*BT[4,2] + k[2]*BT[4,3] + k[3]*BT[4,4])
    elif j==4:
        ys = h*(k[0]*BT[5,1] + k[1]*BT[5,2] + k[2]*BT[5,3] + k[3]*BT[5,4] + k[4]*BT[5,5])
    return ys

@njit
def k_calcs_jit(x,y,h,BT, j, k, xl) :
        
    def dens(x): 
        den = eps * ((1.0 + xl**2)/(1.0 + x**2))
        return den
    
    def func(x,T):
        #den = eps*(1.0 +(rp[i*int(len(rp)/lt)])**2)/(1.0 + (x)**2)  #extra factor of dens withdrawn, j changed to i
        den = dens(x)
        nondim = (r0cgs/RME)*solid_dens 
        c1= 0.153536 #factor that must be an amalgamation of many constants
        B0 = np.log(0.5*(T**2)*(T+2.0)) + (1.0 +(T**2.0)/8.0 - (2.0*T+1)*np.log(2))/((T +1.0)**2.0) #factor from formula
        B = B0 -2.0*np.log(I/(RME*M_fac)) #another factor from formula
        dxde_gas = -nondim*den*c1*zovera*B/(2.0*T) #final stopping power form
        return dxde_gas
    
    if j ==0 :
        k_out = func(x,y)
    elif j==1:
        k_out = func(x + BT[1,0]*h, y - h*(k[0]*BT[1,1])) # should be + not - but need to account for negative dr 
    elif j==2:
        k_out = func(x + BT[2,0]*h, y - h*(k[0]*BT[2,1] + k[1]*BT[2,2]))
    elif j==3:
        k_out = func(x + BT[3,0]*h, y - h*(k[0]*BT[3,1] + k[1]*BT[3,2] + k[2]*BT[3,3]))
    elif j==4:
        k_out = func(x + BT[4,0]*h, y - h*(k[0]*BT[4,1] + k[1]*BT[4,2] + k[2]*BT[4,3] + k[3]*BT[4,4]))
    elif j==5:
        k_out = func(x + BT[5,0]*h, y - h*(k[0]*BT[5,1] + k[1]*BT[5,2] + k[2]*BT[5,3] + k[3]*BT[5,4] + k[4]*BT[5,5])) 
    return k_out 

@njit
def ys_calc_jit(h,k,BT,j):
    if j == 0:
        ys = h*(k[0]*BT[1,1])
    elif j ==1:
        ys = h*(k[0]*BT[2,1] + k[1]*BT[2,2])
    elif j==2:
        ys = h*(k[0]*BT[3,1] + k[1]*BT[3,2] + k[2]*BT[3,3])
    elif j ==3: 
        ys = h*(k[0]*BT[4,1] + k[1]*BT[4,2] + k[2]*BT[4,3] + k[3]*BT[4,4])
    elif j==4:
        ys = h*(k[0]*BT[5,1] + k[1]*BT[5,2] + k[2]*BT[5,3] + k[3]*BT[5,4] + k[4]*BT[5,5])
    return ys
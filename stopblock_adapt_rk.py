import numpy as np 
import scipy.interpolate as spint 
from gen_var import *
from electron import RME as RME_e 
import os 
import rkf45_mod
import getpot as gp

def stop(E, i, particle, savedir,pot):
    from electron import RME, M_fac
    from gen_var import zovera, delta, I, le  
    le_mid = len(E)
    def stop(x,T):
        den = eps*(1.0 +(rp[i*int(len(rp)/lt)])**2)/(1.0 + (rc[i*int(len(rc)/lt)]-x)**2)  #extra factor of dens withdrawn, j changed to i
        nondim = (r0cgs/RME)*solid_dens #non-dimensionalisation factor, in cgs to work with formula,this eliminates units of c1
        c1= 0.153536 #factor that must be an amalgamation of many constants
        B0 = np.log(0.5*(T**2)*(T+2.0)) + (1.0 +(T**2.0)/8.0 - (2.0*T+1)*np.log(2))/((T +1.0)**2.0) #factor from formula
        B = B0 -2.0*np.log(I/(RME*M_fac)) #another factor from formula
        dxde_gas = -nondim*den*c1*zovera*B/(2.0*T) #final stopping power form
        return dxde_gas 
    
    def dens(x,i):
        from gen_var import eps,rp,rc 
        den = eps*((1.0 + rp[i]**2)/(1.0 + (rc[i] - x)**2))
        return den 

    for j in range(le_mid):
        en = [E[j]/(RME*M_fac)] #non-dimensionalising energy for use in STOP function 
        x = 0.0
        x_arr = np.asarry(x)
        en_arr = np.asarray(en)
        x_new, en_new, h, err = rkf45_mod.rkf(x,en,h,stop,rc[i],rc[i] - rp[i])
        if np.isnan(en_new) == 'true':
            pass
        else:
            pot_diff = gp.getpot(pot,h,x_arr[-1])
            en_new += pot_diff
            en_arr = np.append(en_arr, en_new)
            x_arr = np.append(x_arr, x_new)
        en_arr[:] *= RME*M_fac
        en_prof = np.append(x_arr, en_arr)
        en_prof = np.reshape(en_prof, (int(0.5*len(x_arr)), 2),'F')
        np.save(os.path.join(savedir,'EvsR_'+particle +'_E0%d.npy' % (j)), en_prof)
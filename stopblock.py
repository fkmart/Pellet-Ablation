
import numpy as np
import os
from gen_var import *
from electron import RME as RME_e
#from numba import jit, float64
#import stop_func


def stopblock(E,r, i, particle):
    from electron import RME, M_fac
    from gen_var import zovera, delta, I, le  
    le_mid = len(E)
    print('time is', i)
    lr = len(r) #lr will change with time to maintain spatial resolution, must be recalculated at every time
    
    #if (particle=='electron'):
                
    def stop(x,T,j):
        den = eps*(1.0 +(rp[i*int(len(rp)/lt)])**2)/(1.0 + (rc[i*int(len(rc)/lt)]-x)**2)  #extra factor of dens withdrawn, j changed to i
        nondim = (r0cgs/RME)*solid_dens #non-dimensionalisation factor, in cgs to work with formula,this eliminates units of c1
        c1= 0.153536 #factor that must be an amalgamation of many constants
        B0 = np.log(0.5*(T**2)*(T+2.0)) + (1.0 +(T**2.0)/8.0 - (2.0*T+1)*np.log(2))/((T +1.0)**2.0) #factor from formula
        B = B0 -2.0*np.log(I/(RME*M_fac)) #another factor from formula
        dxde_gas = -nondim*den*c1*zovera*B/(2.0*T) #final stopping power form
        return dxde_gas 
        
    for j in range(0, le_mid):
          k = 0
          #r = np.arange(0, rc[j*int(len(rc)/lt)] - rp[j*int(len(rc)/lt)],dr)  #array of zero to pellet surface through cloud
          en = [E[j]/(RME*M_fac)] #non-dimensionalising energy for use in STOP function 
          k1,k2,k3,k4 = 0,0,0,0 #RK4 algorithm terms
          z = np.arange(0,lr-1,1) #helps to determine r index in while loop with y counter
          y = 0 #reset y counter to 0 after every energy
          while en[y] > I/(RME*10**6) and en[y]+0.5*k1 > I/(RME*10**6) and en[y]+0.5*k2 > I/(RME*10**6) and en[y] + k3 >I/(RME*10**6) and y< lr -1:
            energy = en[y] + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)        #calculate new energy        
            k = z[y] #legacy code from first instance
            en.append(energy)            #en is the accumulated energy profile
            k1 = dr*stop(r[k+1],en[y],j) #next four lines are RK4 
            k2 = dr*stop(r[k+1]+dr*0.5, en[y]+k1*0.5,j)
            k3 = dr*stop(r[k+1] + dr*0.5, en[y]+k2*0.5,j)
            k4 = dr*stop(r[k+1] + dr, en[y] + k3,j)              
            y = y+1 # y counter increases
                  
          else:
              en = en[1:] # remove first element of en, dobules up on initial energy
              #print('lr is', lr, 'and the terminal y value is', y)
              for x in range(0, len(en)):
                  en[x] *= RME*M_fac
              #en[:] *= RME*M_fac # converts energy back to eV 
              ener_prof = np.append(en, r[0:k+1]) # add the spacial points to ener_prof as final output
              ener_prof = np.reshape(ener_prof, (int(len(ener_prof)/2),2), 'F') # reshape array for two columns                  
              np.save(os.path.join('raw_multi_static_outputs','EvsR_'+particle +'_t%d_E0%d.npy' % (i,j)), ener_prof) # save files as npy for analysis

 
"""This is the defining stopblock function - as in, block of stopping code
The following function (appended with phi) includes a consideration of the potential field the particles are moving through also.
Neutral density is defined across the relevant points on the general grid, the stoping functino for use with RK4 is included too
THe RK4 solver is defined here explicitly because the sheer number of possibilities for failure with log functions etc made this easier
"""

def stopblock_phi(E,r, i, den,phi, savedir):
        from electron import RME, M_fac
        from gen_var import zovera, delta, I, le  
        le_mid = len(E)
        particle = 'electron'
        lr = len(r) #lr will change with time to maintain spatial resolution, must be recalculated at every time

        def dens(x):
            if x < rp[i]:
                dens = 0.0
            elif x > rc[i]:
                dens = 0.0
            else:
                dens = 0.01 * ((1.0 + rp[i]**2)/(1.0 + x**2))
            return dens
        def bethe(x,T):
            #den = eps*(1.0 +(rp[i*int(len(rp)/lt)])**2)/(1.0 + (x)**2)  #extra factor of dens withdrawn, j changed to i
            den = dens(x)
            nondim = (r0cgs/RME)*solid_dens 
            c1= 0.153536 #factor that must be an amalgamation of many constants
            B0 = np.log(0.5*(T**2)*(T+2.0)) + (1.0 +(T**2.0)/8.0 - (2.0*T+1)*np.log(2))/((T +1.0)**2.0) #factor from formula
            B = B0 -2.0*np.log(I/(RME*M_fac)) #another factor from formula
            dxde_gas = -nondim*den*c1*zovera*B/(2.0*T) #final stopping power form
            return dxde_gas 
        
        for j in range(0, le_mid):
              en = [E[j]/(RME*M_fac)] #non-dimensionalising energy for use in STOP function 
              k1,k2,k3,k4 = 0,0,0,0 #RK4 algorithm terms
              #z = np.arange(0,lr-1,1) #helps to determine r index in while loop with y counter
              y = 0 #reset y counter to 0 after every energy
              dphi = 0.0
              while en[y] > I/(RME*10**6) and en[y]+0.5*k1 > I/(RME*10**6) and en[y]+0.5*k2 > I/(RME*10**6) and en[y] + k3 >I/(RME*10**6) and y< lr -1 and en[y] + dphi > I/(RME*M_fac):
                if (y <lr-2):
                    phi1 = phi[lr-y -1]
                    phi2 = phi[lr - y - 2]
                    dphi = phi2 - phi1
                else:
                    dphi = 0.0
                energy = en[y] + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4) #calculate new energy        
                energy += dphi
                en.append(energy)            #en is the accumulated energy profile
                k1 = dr*bethe(r[-y-1],en[y]) #next four lines are RK4 
                k2 = dr*bethe(r[-y-1]+dr*0.5, en[y]+k1*0.5)
                k3 = dr*bethe(r[-y-1] + dr*0.5, en[y]+k2*0.5)
                k4 = dr*bethe(r[-y-1] + dr, en[y] + k3)              
                y = y+1 # y counter increases
                  
              else:
                  en = en[1:] # remove first element of en, dobules up on initial energy
                  for x in range(0, len(en)):
                      en[x] *= RME*M_fac
                  ener_prof = np.append(en, r[0:y]) # add the spacial points to ener_prof as final output
                  ener_prof = np.reshape(ener_prof, (int(len(ener_prof)/2),2), 'F') # reshape array for two columns                  
                  np.save(os.path.join(savedir,'EvsR_'+particle +'_E0%d.npy' % (j)), ener_prof) # save files as npy for analysis
            
"""
    elif(particle == 'proton'):
        
        def stop(x, T, j):
            den = eps*(1.0 +(rp[i*int(len(rp)/lt)])**2)/(1.0 + (rc[i*int(len(rc)/lt)]-x)**2)  #extra factor of dens withdrawn, j changed to i
            nondim = (r0cgs/RME)*solid_dens #non-dimensionalisation factor, in cgs to work with formula
           
            c1= 0.153536 #factor that must be an amalgamation of many constants
            
            #B0 = np.log(0.5*(T**2)*(T+2.0)) + (1.0 +(T**2.0)/8.0 - (2.0*T+1)*np.log(2))/((T +1.0)**2.0) #factor from formula
            B = -2.0*np.log(I/(RME*M_fac)) #another factor from formula
            numer = 2.0*T*RME_e/RME
            denom = I/(RME*M_fac)
            logarg = numer/denom
            B = 2.0*np.log(logarg)
            dxde_gas = -nondim*den*c1*zovera*B/(2.0*T) #final stopping power form
            
            return dxde_gas 
        
        for j in range(0, le_mid):
              #r = np.arange(0, rc[j*int(len(rc)/lt)] - rp[j*int(len(rc)/lt)],dr)  #array of zero to pellet surface through cloud
              en = [E[j]/(RME*M_fac)] #non-dimensionalising energy for use in STOP function
              
              
              k1 = dr*stop(r[0],en[0],j) #next four lines are RK4 
              k2 = dr*stop(r[0]+dr*0.5, en[0]+k1*0.5,j)
              k3 = dr*stop(r[0] + dr*0.5, en[0]+k2*0.5,j)
              k4 = dr*stop(r[0] + dr, en[0] + k3,j)    
              #k1,k2,k3,k4 = 0,0,0,0 #RK4 algorithm terms
              z = np.arange(0,lr-1,1) #helps to determine r index in while loop with y counter
              y = 0 #reset y counter to 0 after every energy
              while en[y] > I/(RME*10**6) and en[y]+0.5*k1 > I/(RME*10**6) and en[y]+0.5*k2 > I/(RME*10**6) and en[y] + k3 >I/(RME*10**6) and y< lr -1 and 2.0*en[y]*RME_e/RME>I/(RME*M_fac):
                  energy = en[y] + (1.0/6.0)*(k1 + 2*k2 + 2*k3 + k4)        #calculate new energy         
                  k = z[y] #legacy code from first instance
                  en.append(energy)            #en is the accumulated energy profile
                  k1 = dr*stop(r[y+1],en[y],j) #next four lines are RK4 
                  k2 = dr*stop(r[y+1]+dr*0.5, en[y]+k1*0.5,j)
                  k3 = dr*stop(r[y+1] + dr*0.5, en[y]+k2*0.5,j)
                  k4 = dr*stop(r[y+1] + dr, en[y] + k3,j)              
                  y = y+1 # y counter increases
                  
              else:
                  #en = en[1:] # remove first element of en, dobules up on initial energy
                  #print(en)
                  #print('lr is', lr, 'and the terminal y value is', y)
                  for x in range(0, len(en)): 
                      en[x] = en[x]*RME*M_fac # loop through en array amd convert to eV
                  #print (en)
                  #print('JUst energy array has shape', np.shape(en))   
                  ener_prof = np.append(en, r[0:y+1]) # add the spacial points to ener_prof as final output
                  #print(ener_prof)
                  #print('Energy profile with radial points is', np.shape(ener_prof))
                  ener_prof = np.reshape(ener_prof, (int(len(ener_prof)/2),2), 'F') # reshape array for two columns 
                  
                  np.save(os.path.join('Outputs_phi','EvsR_'+particle +'_t%d_E0%d.npy' % (i,j)), ener_prof) # save files as npy for analysis
         #return en"""
         #np.save('OUT_CALC_IMP_elecstop_fs_1keVbroad_t%d_E0_%d.npy' % (j,i), en)
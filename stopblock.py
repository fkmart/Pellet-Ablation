
import numpy as np
import os
from gen_var import *
from electron import RME as RME_e
import dphi_calc
from gen_var import trunc_fac
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


def stopblock_phi_mod(E,r, i, den,phi, savedir):
        from electron import RME, M_fac
        from gen_var import zovera, delta, I, le, dr
        import rk45_mod as rk 
        le_mid = len(E)
        particle = 'electron'
        lr = len(r) #lr will change with time to maintain spatial resolution, must be recalculated at every time
        dr = r[1] - r[0]
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
            thing = 2.0*np.log(I/(RME*M_fac))
            B = B0 -2.0*np.log(I/(RME*M_fac)) #another factor from formula
            dxde_gas = -nondim*den*c1*zovera*B/(2.0*T) #final stopping power form
            return dxde_gas 
        
        for j in range(0, le_mid):
              print(j)
              if j ==100:
                  print('here')
              en = [E[j]/(RME*M_fac)] #non-dimensionalising energy for use in STOP function 
              k1,k2,k3,k4 = 0,0,0,0 #RK4 algorithm terms
              #z = np.arange(0,lr-1,1) #helps to determine r index in while loop with y counter
              y = 0 #reset y counter to 0 after every energy
              dphi = 0.0
              x = rc[i]
              x_tot = np.asarray(x)
              I_non = I/(RME*M_fac)
              energy = en[0]
              dr = 1e-3
              while en[y] > I_non and en[y]+0.5*k1 > I_non and en[y]+0.5*k2 > I_non and en[y] + k3 >I_non and x>rp[i] and en[y] + dphi > I_non:
                if x > rp[i]:#if (y <lr-2):
                    #phi1 = phi[-y ]
                    #phi2 = phi[ - y - 1]
                    #dphi = phi2 - phi1
                    dphi = dphi_calc.dphi_calc(x+np.abs(dr),x,r_grid, phi)
                else:
                    dphi = 0.0
                #energy += dphi
                            #en is the accumulated energy profile   
                 
                
                energy, k1, k2, k3, k4 = rk.rk(x, en[-1],bethe,-dr) 
                x -= dr
                energy += dphi
                en.append(energy)
                x_tot = np.append(x_tot,x)
                y = y+1 # y counter increases
                if x - np.abs(dr) < rp[i]:
                    dr = x - rp[i]
                else:
                    pass
              else:
                  en = en[:] # remove last element of en #THIS IS THE POSSIBLE PROBLEM
                  x_tot = x_tot[:]
                  if j ==161:
                      print('here')
                  for s in range(0, len(en)):
                      en[s] *= RME*M_fac
                  ener_prof = np.append(en, x_tot) # add the spacial points to ener_prof as final output
                  ener_prof = np.reshape(ener_prof, (int(len(ener_prof)/2),2), 'F') # reshape array for two columns                  
                  np.save(os.path.join(savedir,'EvsR_'+particle +'_E0%d.npy' % (j)), ener_prof) # save files as npy for analysis           

def stopblock_phi_mod_rkf(E,r_grid, i,phi, savedir):
    from electron import RME, M_fac
    from gen_var import zovera, delta, I, le, eps, rp, r0cgs, solid_dens,rc, style
    import rkf45_CSDA as rkf 
    import step_calc as sc

    if style =='many':
        from gen_var import rp_hr as rp
        from gen_var import rc_hr as rc
    else:
        pass

    le_mid = len(E)

    particle = 'electron'

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

    def dens(x): 
        den = eps * ((1.0 + rp[i]**2)/(1.0 + x**2))
        return den

    def bethe(x,T):
        #den = eps*(1.0 +(rp[i*int(len(rp)/lt)])**2)/(1.0 + (x)**2)  #extra factor of dens withdrawn, j changed to i
        den = dens(x)
        nondim = (r0cgs/RME)*solid_dens 
        c1= 0.153536 #factor that must be an amalgamation of many constants
        B0 = np.log(0.5*(T**2)*(T+2.0)) + (1.0 +(T**2.0)/8.0 - (2.0*T+1)*np.log(2))/((T +1.0)**2.0) #factor from formula
        thing = 2.0*np.log(I/(RME*M_fac))
        B = B0 -2.0*np.log(I/(RME*M_fac)) #another factor from formula
        dxde_gas = -nondim*den*c1*zovera*B/(2.0*T) #final stopping power form
        return dxde_gas

    for j in range(0, le_mid):
        en = [E[j]/(RME*M_fac)] #non-dimensionalising energy for use in STOP function 
        k1,k2,k3,k4,k5,k6 = 0,0,0,0,0,0 #RKF4 algorithm terms
        #z = np.arange(0,lr-1,1) #helps to determine r index in while loop with y counter
        y = 0 #reset y counter to 0 after every energy
        dphi = 0.0
        tot_pot = 0.0
        x = rc[i]
        x_tot = np.asarray(x)
        I_non = I/(RME*M_fac)
        y_steps = np.zeros(5)
        dr = 1e-3
        while en[y] > trunc_fac*I_non and en[y] - y_steps[0] >trunc_fac*I_non and en[y] -y_steps[1] > trunc_fac*I_non and en[y] - y_steps[2] > trunc_fac*I_non and en[y] - y_steps[3] > trunc_fac*I_non and en[y] - y_steps[4] > trunc_fac*I_non and x>rp[i] and en[y] + dphi > trunc_fac*I_non:
          #en[y]-dr*(k1*BT[1,1]) > I_non and en[y]-dr*(k1*BT[2,1] + 
          #k2*BT[2,2]) > I_non and en[y] - dr*(k1*BT[3,1] + k2*BT[3,2] + 
          #k3*BT[3,3]) >I_non and en[y] - dr*(k1*BT[4,1] + 
          #k2*BT[4,2] + k3*BT[4,3] + k4*BT[4,4]) > I_non and en[y] - dr*(k1*BT[5,1] + k2*BT[5,2] + 
          #k3*BT[5,3] + k4*BT[5,4] + k5*BT[5,5]) > I_non and x>rp[i] and en[y] + dphi > I_non:
        
            x, energy,dr, err, new_k = rkf.rkf(x, en[-1],dr, bethe, BT, BTS, rc[i], rp[i],I_non)
            k1,k2,k3,k4,k5,k6 = new_k[0], new_k[1], new_k[2], new_k[3], new_k[4], new_k[5]
            dphi = dphi_calc.dphi_calc(x+np.abs(dr),x,r_grid, phi)
            tot_pot += dphi
            energy += dphi
            en.append(energy)            #en is the accumulated energy profile   
            x_tot = np.append(x_tot,x)
            y = y+1 # y counter increases
            k = [k1,k2,k3,k4,k5,k6]
            for p in range(0,5):
                y_steps[p] = sc.ys_calc(dr,k,BT,p)
        else:
              en = en[:] 
              x_tot = x_tot[:]
              en_check = np.asarray(en[-200:]) * RME*M_fac
              tot_pot *= RME*M_fac
              x_check = x_tot[-200:]
              for s in range(0, len(en)):
                  en[s] *= RME*M_fac
              ener_prof = np.append(en, x_tot) # add the spacial points to ener_prof as final output
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
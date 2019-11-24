import numpy as np
from electron import e_bar, e_dist, ener_res, dist_calc 
from gen_var import rp, rc , t_start, sig, pel_pot , lp, I 
import gauss_test_pot as gtp
import bethe_analysis_functions as baf 
import os 

direc = os.getcwd() + '/pot_stop'
subdir = '/data/'

filelist_raw = [ f for f in os.listdir(direc + subdir ) if f.endswith(".npy") ] # deletes the .npy files of raw data
for f in filelist_raw:
    os.remove(os.path.join(direc + subdir, f))

#No need for grid, but why not 
t = t_start
r1 = rp[t]
r2 = rc[t]
r_int = np.linspace(r1, r2, 602) # 602 selected to give same resolution as bethe simulations
lr = len(r_int)
mid = int(len(r_int)*0.5)

e_mid, e_bins, mb = dist_calc(e_dist,ener_res,e_bar)

for a in range(0, len(sig)):
    print('The ' + str(a) + 'th iteration of the standard deviation')
    for p in range(1,lp):
        pot = gtp.gauss_func(pel_pot[p],sig[a],r_int[mid], r_int)
        for i in range(0, len(e_mid)):
            en = e_mid[i]
            energy = []
            energy.append(en)
            y = 0
            while (y < lr -1) and (en > I):
                pot1 = pot[lr - 1 -y]
                pot2 = pot[lr - 2 - y]
                dphi = pot2 - pot1
                y +=1
                en += dphi
                energy.append(en)
            x_stop = r_int[-y] 
            traj = np.flip(r_int[-1-y:], axis = 0)
            ener_prof = energy 
            ener_prof = np.append(energy, traj)
            ener_prof = np.reshape(ener_prof, (int(len(ener_prof)*0.5),2),'F')
            np.save(direc + subdir + 'prof%d_pot%d_sig%.2f'%(i,int(pel_pot[p]),sig[a]), ener_prof)
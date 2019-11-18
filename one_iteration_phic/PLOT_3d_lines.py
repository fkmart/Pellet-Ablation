import numpy as np 
import matplotlib.pyplot as plt 
from mpl_toolkits import mplot3d
from gen_var import p_inc, pel_pot
import os 
direc = os.getcwd() 

load_dir = direc + '/one_iteration_phic/'

k = [0.5,1.0,2.0] 
for a in range(0,len(k)):
    for p in range(0, len(p), p_inc):

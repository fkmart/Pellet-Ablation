# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 16:42:23 2020

@author: Kyle
"""

import os 
import shutil

direc = os.getcwd() 

newpath = os.path.join(direc, 'many_iteration','neutral','500eV','analysed_outputs')

d = [x[0] for x in os.walk(newpath)][1:]

shutil.rmtree(direc + os.sep + 'test1')
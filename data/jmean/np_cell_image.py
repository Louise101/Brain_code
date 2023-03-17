#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 15:57:55 2022

@author: lf58
"""

import matplotlib.pyplot as plt
import numpy as np
import math as math
import matplotlib.image as mpimg
import csv



depth=[]
n = 100
fd = open('np_cell_full.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.int32).reshape(n,n,101)
fd.close()
data=np.flipud(datain)

zlength=2 #cm
znumber=200
zvoxel_length=zlength/znumber #cm

dslice=data[:n,0:n,50] #use data for voxel number axis and datain for axis produced in python like depth. 


#for i in range(znumber):
#	idepth=(i*zvoxel_length) #cm
#	depth.append(idepth)
	

#image=plt.imshow(dslice,extent=[min(depth),max(depth),min(depth),max(depth)])
plt.imshow(dslice)
plt.xlabel('x (cm)') 
plt.ylabel('z (cm)') 
plt.title('skin layers')
plt.colorbar()
plt.show()

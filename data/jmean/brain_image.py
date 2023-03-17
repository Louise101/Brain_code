#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar  4 16:14:03 2022

@author: lf58
"""


import matplotlib.pyplot as plt
import numpy as np
import math as math
import matplotlib.image as mpimg
import csv



depth=[]

fd = open('rhokap.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(150,206,236)
fd.close()
data=np.flipud(datain)

#zlength=2 #cm
#znumber=200
#zvoxel_length=zlength/znumber #cm

#dslice=datain[13:162,50:255,132] #use data for voxel number axis and datain for axis produced in python like depth. 
#dslice=datain[92,50:255,20:255] #use data for voxel number axis and datain for axis produced in python like depth. 

dslice=datain[75,:206,:236]

#APPLY TO MATLAB GRIDS AND INTO FORTRAN - also try transpose . ..


#for i in range(znumber):
#	idepth=(i*zvoxel_length) #cm
#	depth.append(idepth)
	

#image=plt.imshow(dslice,extent=[min(depth),max(depth),min(depth),max(depth)])
plt.imshow(dslice)
plt.xlabel('x (cm)') 
plt.ylabel('z (cm)') 
plt.title('Greymatter')
plt.colorbar()
plt.show()
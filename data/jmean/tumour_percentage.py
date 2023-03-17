#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep  5 11:28:13 2022

@author: lf58
"""

import matplotlib.pyplot as plt
import numpy as np
import math as math
import matplotlib.image as mpimg
import csv
import matplotlib.colors as color
from matplotlib import pyplot as plt, cm

#[z,y,x]
depth=[]

#fd = open('1rhokap.dat', 'rb')
#fd = open('1tum_kill_code12.dat', 'rb')
#datain=np.fromfile(file=fd, dtype=np.double).reshape(155,231,250)
#fd.close()
#data_1=datain#np.flip(datain)

#fd = open('60tum_kill_code12.dat', 'rb')
#datain=np.fromfile(file=fd, dtype=np.double).reshape(155,231,250)
#fd.close()
#data=datain#np.flip(datain)


fd = open('tumkill_slice16.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(1200,231,250)
fd.close()
tum=datain#np.flip(datain)

data_1=tum[0,:,:]
data=tum[360,:,:]

init_tum=0.
final_tum=0.
#get totals of tumour grids
for i in range (len(data)):
    for j in range (len(data[0])):
        #for k in range(len(data[0][0])):
            init_tum=init_tum + data_1[i][j]#[k]
            final_tum=final_tum + data[i][j]#[k]

percentage_left=(final_tum/init_tum)*100. 
percentage_killed=100-percentage_left

print(init_tum, final_tum, percentage_left, percentage_killed)
            
            


dslice_1=data_1#[77,:,:]
dslice=data#[77,:,:]


for i in range (len(dslice)):
    for j in range (len(dslice[0])):
            if (dslice_1[i][j]-dslice[i][j] != 0.):
                dslice[i][j]= np.nan
            if (dslice_1[i][j]== 10):
                dslice[i][j]=0.
            
                



#150 - (75-13)=88
#206-(200-50)=56
#236-(220-20)=36

#APPLY TO MATLAB GRIDS AND INTO FORTRAN - also try transpose . ..

ylength=5.38 #cm
ynumber=231
yvoxel_length=ylength/ynumber #cm

ydepth=[]

for i in range(ynumber):
	idepth=(i*yvoxel_length) #cm
	ydepth.append(idepth)
	
xlength=5.82 #cm
xnumber=250
xvoxel_length=xlength/xnumber #cm

xdepth=[]

for i in range(xnumber):
	idepth=(i*xvoxel_length) #cm
	xdepth.append(idepth)


cmap=cm.get_cmap()
cmap.set_bad('red') 
    
image=plt.imshow(dslice,extent=[min(xdepth),max(xdepth),min(ydepth),max(ydepth)])#,norm=color.LogNorm())	
#plt.imshow(dslice)
#plt.xlabel('x (cm)') 
#plt.ylabel('y (cm)') 
#plt.title('1 Sec - Heterogeneous')
plt.colorbar()
plt.savefig('Tumour_600_thresh560.png',format='png',dpi=1200, bbox_inches='tight')
plt.show()
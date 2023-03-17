#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 11:32:07 2023

@author: lf58
"""

#need to renormalise after interploation and ballon max = 100 not 1!!

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as color
from scipy.ndimage import zoom

fd = open('tumour.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(150,206,236)
fd.close()
tumour=np.flip(datain)

print(tumour.max())



tumour_zoom=zoom(tumour,3, mode='nearest')

#multiplying by 4 - so voxel edge size is now 0.7mm/4.

print(len(tumour_zoom))

#dslice_cut=tumour_zoom[273,35:251,65:301]
#dslice_cut=tumour_zoom[190:350,108,65:301]


#print(dslice_cut[100,:])

#print(tumour[1,1,1])
maxi=[]
for i in range (len(tumour_zoom)):
    for j in range (len(tumour_zoom[0])):
        for k in range(len(tumour_zoom[0][0])):
            if(tumour_zoom[i,j,k] < 0.1):
                tumour_zoom[i,j,k]=0.
                
      #      if(tumour_zoom[i,j,k] > 1.):
      #          tumour_zoom[i,j,k]=1.
                
            
                
#for i in range (len(tumour_zoom)):
 #   for j in range (len(tumour_zoom[0])):
  #      for k in range(len(tumour_zoom[0][0])):
   #         if(tumour_zoom[i,j,k] > 90.):
    #            tumour_zoom[i,j,k]=100.
     #       else:
      #          tumour_zoom[i,j,k]=0.
          
                
                
            #if(tumour_zoom[i,j,k] <100.):
             #   tumour_zoom[i,j,k]=0.


dslice=tumour_zoom[273,:,:]

dslice_cut=tumour_zoom[273,35:251,65:301]
#dslice_cut=tumour_zoom[190:350,108,65:301]


print(dslice_cut[100,:])
dslice_org=tumour[80,:83,20:100]

#dslice=tumour[65:115,25,20:100]

#tum_cut=tumour[65:115,:83,20:100]

# dimensions = 50,83,80 

#dslice_cut=tum_cut[:,:,90]

tumour_zoom_cut=tumour_zoom[195:350,25:256,55:305]

dtzc=tumour_zoom_cut[152,:,:]

plt.imshow(dslice_cut)
#plt.xlabel('x (mm)') 
#plt.ylabel('y (mm)') 
#plt.title('Jmean')
plt.colorbar()
#plt.savefig('cav_flu.png',format='png',dpi=1200, bbox_inches='tight')
plt.show()

flat=tumour_zoom_cut.flatten()
    
f=np.savetxt('tum_hires2.txt', flat)
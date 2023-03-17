#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon 9th Jan 2023

@author: lf58
"""

import matplotlib.pyplot as plt
import matplotlib.colors as color
import numpy as np
import math as math
import matplotlib.image as mpimg
import csv

#[z,y,x]
depth=[]

fd = open('jmean_brain_code26.dat', 'rb')
#fd = open('60o21_g21_ox39_thresh560.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
fd.close()
data=datain#np.flip(datain)

#fd = open('600tum_kill_g21_ox39_thresh560.dat', 'rb')
#datain=np.fromfile(file=fd, dtype=np.double).reshape(150,206,236)
#fd.close()
#data_tum=np.flip(datain)


fd = open('rhokap_brain_code26.dat', 'rb')
datain_rho=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
fd.close()
rho=datain_rho#np.flip(datain_rho)

for i in range (len(rho)):
  for j in range (len(rho[0])):
      for k in range(len(rho[0][0])):
         if (rho[i][j][k]==10.001 or rho[i][j][k]==0.020001 or rho[i][j][k]==0.39 or rho[i][j][k]==0. or rho[i][j][k]==0.003):
            data[i][j][k]=0.
            rho[i][j][k]=0.

#for i in range (len(data)):
 #   for j in range (len(data[0])):
  #      for k in range(len(data[0][0])):
   #         if (rho[i][j][k]==11.0019 or rho[i][j][k]==0.):
    #           data[i][j][k]=0.
    
    
#ax=plt.figure().add_subplot(projection='3d')
#ax.voxels(rho)
#plt.show()
                
#data[30,40,40]=0.vpn.st-andrews.ac.uk


#zlength=2 #cm
#znumber=200
#zvoxel_length=zlength/znumber #cm

#starting centre of tumour at 75,200,220 in matlab original data cube of 176,320,320 . . .
#in the cutversion this is: 88,56,36

#dslice=datain[13:162,50:255,132] #use data for voxel number axis and datain for axis produced in python like depth. 
#dslice=datain[92,50:255,20:255] #use data for voxel number axis and datain for axis produced in python like depth. 

#Slices where pure tumour is seen (when rhokap set to 1,2,3 for each tissue type)
#so possibly assume these points are at the centre of the cavity? 
#dslice=data[91,:206,:236]
#dslice=data[:150,35,:236]
#dslice=data[:150,:206,58]

print(np.max(data))
#print(data[0][0])

#FIND THE 500 USIG ELEMENT SEARCH  . . .

#dslice=data[92,10:50,30:85]
#dslice=data[91,:65,:110]# + rho[91,:206,:236]
dslice=rho[77,:,:]# +rho[91,:,:]
#dslice=data[91,:75,:100]
#dslice=rho[75,:206,:236]
#dslice_cut=dslice[20:80,10:60]
#print(dslice_cut[50])

#block=rho[77,128:173,174:219]

#dslice=block



#dslice=data[91,:,:]
#dslice=data[:,35,:]
# dslice=data[:,:,58] 

#dslice=data[91,:,:]



with open('brain__rho_jmean_slice.csv','w+')as my_csv:
    csvWriter=csv.writer(my_csv,delimiter=',')
    csvWriter.writerows(dslice)

#print(list(dslice))

#150 - (75-13)=88
#206-(200-50)=56
#236-(220-20)=36

#APPLY TO MATLAB GRIDS AND INTO FORTRAN - also try transpose . ..


#for i in range(znumber):
#	idepth=(i*zvoxel_length) #cm
#	depth.append(idepth)

ylength=14.42 #cm
ynumber=206
yvoxel_length=ylength/ynumber #cm

ydepth=[]

for i in range(ynumber):
	idepth=(i*yvoxel_length) #cm
	ydepth.append(idepth)
	
xlength=16.52 #cm
xnumber=236
xvoxel_length=xlength/xnumber #cm

xdepth=[]

for i in range(xnumber):
	idepth=(i*xvoxel_length) #cm
	xdepth.append(idepth)
    
#image=plt.imshow(dslice,extent=[min(xdepth),max(xdepth),min(ydepth),max(ydepth)])#,norm=color.LogNorm())	

#image=plt.imshow(dslice,extent=[min(depth),max(depth),min(depth),max(depth)])
plt.imshow(dslice)#, norm=color.LogNorm())
plt.xlabel('x (mm)') 
plt.ylabel('y (mm)') 
#plt.title('Jmean')
plt.colorbar()
plt.savefig('cav_flu.png',format='png',dpi=1200, bbox_inches='tight')
plt.show()
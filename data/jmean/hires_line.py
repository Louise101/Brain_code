#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 26 17:05:16 2023

@author: lf58
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as color

fd = open('rhokap_brain_code26.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
fd.close()
rhokap=datain#np.flip(datain)

fd = open('albedo_brain_code26.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
fd.close()
albedo=datain#np.flip(datain)

#fd = open('hgg_brain_code25.dat', 'rb')
#datain=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
#fd.close()
#hgg=datain#np.flip(datain)

rhoslice=rhokap[77,:,:]

#fd = open('jmean_brain_g21_ox39_hires_ng.dat', 'rb')
#fd = open('600o23_g21_ox39_thresh560_con023_hires.dat', 'rb')
#fd = open('600o21_code7.dat', 'rb')
fd = open('rhokap_brain_code26.dat', 'rb')
fd = open('jmean_brain_code26.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
fd.close()
jmean=datain#np.flip(datain)

for i in range (len(rhokap)):
  for j in range (len(rhokap[0])):
      for k in range(len(rhokap[0][0])):
          if ( rhokap[i][j][k]==10.001 or rhokap[i][j][k]==0.006  ):
            jmean[i][j][k]=0.

depth_x=[]
jmean_x=[]
dose_x=[]
rho_x=[]
alb_x=[]
hgg_x=[]
             
depth_z=[]
jmean_z=[]
dose_z=[]

depth_y=[]
jmean_y=[]
dose_y=[]

for k in range(len(jmean)):
    depth_z.append(k*0.233)
    jmean_z.append(jmean[k][115][125])
    dose_z.append(jmean[k][115][125])#*(558/1000)) #J/cm2

for j in range(len(jmean[0])):
    depth_y.append(j*0.233)
    jmean_y.append(jmean[77][j][125])
    dose_y.append(jmean[77][j][125] )#*(558/1000))
    
for i in range(len(jmean[0][0])):
    depth_x.append(i*0.233)
    jmean_x.append((jmean[77][155][i]))
    dose_x.append((jmean[77][155][i])*(558/1000))
    rho_x.append((rhokap[77][155][i]))
    alb_x.append((albedo[77][155][i]))
   # hgg_x.append((hgg[77][153][i]))
    
print(len(rho_x[195:]))
print(rho_x[195:])
print(alb_x[195:])
#print(hgg_x[201:])
    
dslice_xy=jmean[77,:,:]
image=plt.imshow(dslice_xy,extent=[min(depth_x),max(depth_x),min(depth_y),max(depth_y)])#,cmap='Reds',norm=color.LogNorm())	
plt.colorbar()
plt.show()

    
    
#plt.plot(depth_x[:36],dose_x[:36], linewidth='3', label='left')
plt.plot(depth_x[100:],jmean_x[100:], linewidth='3', label='right')

#plt.plot(depth_x[85:105],jmean_x_depth, linestyle='dashed', color='black')
#plt.plot(jmean_x_flu,jmean_x[85:105])
#plt.plot(jmean_x_25,jmean_x[85:105])
#plt.plot(jmean_x_wall,jmean_x[84:105])
#plt.yscale('log')
plt.xlabel('x-depth (mm)') 
plt.ylabel('Fluence Rate (mW/cm2)')  
plt.title('x depth fluence rate -left')
#plt.legend()
#plt.xticks(fontsize=15)
#plt.yticks(fontsize=15)
plt.savefig('flu_wall_x_code1_right.png',format='png',dpi=1200, bbox_inches='tight')
plt.show()
    
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb  8 14:51:11 2023

@author: lf58
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as color

fd = open('rhokap_brain_code14.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
fd.close()
rhokap=datain#np.flip(datain)


fd = open('1o21_code14.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
fd.close()
s01_13=datain#np.flip(datain)

fd = open('30o21_code14.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
fd.close()
s030_13=datain#np.flip(datain)

fd = open('60o21_code14.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
fd.close()
s060_13=datain#np.flip(datain)

fd = open('120o21_code14.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
fd.close()
s0120_13=datain#np.flip(datain)

fd = open('180o21_code14.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
fd.close()
s0180_13=datain#np.flip(datain)

fd = open('240o21_code14.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
fd.close()
s0240_13=datain#np.flip(datain)

fd = open('300o21_code14.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
fd.close()
s0300_13=datain#np.flip(datain)

fd = open('420o21_code14.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
fd.close()
s0420_13=datain#np.flip(datain)

fd = open('480o21_code14.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
fd.close()
s0480_13=datain#np.flip(datain)

#fd = open('540o21_code13.dat', 'rb')
#datain=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
#fd.close()
#s0540_13=datain#np.flip(datain)

#fd = open('600o21_code13.dat', 'rb')
#datain=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
#fd.close()
#s0600_13=datain#np.flip(datain)




#for i in range (len(rhokap)):
 # for j in range (len(rhokap[0])):
  #    for k in range(len(rhokap[0][0])):
   #      if (rhokap[i][j][k]==10.001 or rhokap[i][j][k]==0.020001 or rhokap[i][j][k]==0.39 or rhokap[i][j][k]==0.):
    #        s01_13[i][j][k]=0.

depth_x=[]
#jmean_x=[]
dose_x_1=[]
#rho_x=[]
dose_x_30=[]
dose_x_60=[]
dose_x_120=[]
dose_x_180=[]
dose_x_240=[]
dose_x_300=[]
dose_x_420=[]
dose_x_480=[]
#dose_x_540=[]
#dose_x_600=[]
             
#depth_z=[]
#jmean_z=[]
#dose_z=[]

#depth_y=[]
#jmean_y=[]
#dose_y=[]

#for k in range(len(jmean)):
 #   depth_z.append(k*0.233)
  #  jmean_z.append(jmean[k][115][125])
   # dose_z.append(jmean[k][115][125])#*(558/1000)) #J/cm2

#for j in range(len(s01_13[0])):
    #depth_y.append(j*0.233)
    #jmean_y.append(s01_13[77][j][125])
    #dose_y.append(s01_13[77][j][125] )#*(558/1000))
s0=[]   
for i in range(len(s01_13[0][0])):
    depth_x.append(i*0.233)
   # jmean_x.append((s01_13[77][155][i]))
    dose_x_1.append((s01_13[77][155][i]))#*(558/1000))
    dose_x_30.append((s030_13[77][155][i]))#*(558/1000))
    dose_x_60.append((s060_13[77][155][i]))#*(558/1000))
    dose_x_120.append((s0120_13[77][155][i]))#*(558/1000))
    dose_x_180.append((s0180_13[77][155][i]))#*(558/1000))
    dose_x_240.append((s0240_13[77][155][i]))#*(558/1000))
    dose_x_300.append((s0300_13[77][155][i]))#*(558/1000))
    dose_x_420.append((s0420_13[77][155][i]))#*(558/1000))
    dose_x_480.append((s0480_13[77][155][i]))#*(558/1000))
   # dose_x_540.append((s0540_13[77][155][i]))#*(558/1000))
   # dose_x_600.append((s0600_13[77][155][i]))#*(558/1000))



s0.append((s01_13[77][155][115]))
s0.append((s030_13[77][155][115]))
s0.append((s060_13[77][155][115]))
s0.append((s0120_13[77][155][115]))
s0.append((s0180_13[77][155][115]))
s0.append((s0240_13[77][155][115]))
s0.append((s0300_13[77][155][115]))
s0.append((s0420_13[77][155][115]))
s0.append((s0480_13[77][155][115]))
#s0.append((s0540_13[77][155][115]))
#s0.append((s0600_13[77][155][115]))
    
#dslice_xy=s01_13[77,:,:]
#image=plt.imshow(dslice_xy,extent=[min(depth_x),max(depth_x),min(depth_y),max(depth_y)])#,cmap='Reds',norm=color.LogNorm())	
#plt.colorbar()
#plt.show()
time=[1,30, 60,120,180,240,300,420,480]



    
    
#plt.plot(depth_x[:36],dose_x[:36], linewidth='3', label='left')
plt.plot(depth_x[70:160],dose_x_1[70:160], linewidth='3', label='1s')
plt.plot(depth_x[70:160],dose_x_30[70:160], linewidth='3', label='30s')
plt.plot(depth_x[70:160],dose_x_60[70:160], linewidth='3', label='60s')
plt.plot(depth_x[70:160],dose_x_120[70:160], linewidth='3', label='120s')
plt.plot(depth_x[70:160],dose_x_180[70:160], linewidth='3', label='180s')
plt.plot(depth_x[70:160],dose_x_240[70:160], linewidth='3', label='240s')
plt.plot(depth_x[70:160],dose_x_300[70:160], linewidth='3', label='300s')
plt.plot(depth_x[70:160],dose_x_420[70:160], linewidth='3', label='420s')
plt.plot(depth_x[70:160],dose_x_480[70:160], linewidth='3', label='480s')
#plt.plot(depth_x[70:160],dose_x_540[70:160], linewidth='3', label='540s')
#plt.plot(depth_x[70:160],dose_x_600[70:160], linewidth='3', label='600s')
#plt.plot(depth_x[85:105],jmean_x_depth, linestyle='dashed', color='black')
#plt.plot(jmean_x_flu,jmean_x[85:105])
#plt.plot(jmean_x_25,jmean_x[85:105])
#plt.plot(jmean_x_wall,jmean_x[84:105])
#plt.yscale('log')
plt.xlabel('x-depth (mm)') 
plt.ylabel('Concentration (uM)')  
plt.title('o21')
plt.legend()
#plt.xticks(fontsize=15)
#plt.yticks(fontsize=15)
plt.savefig('o211noox.png',format='png',dpi=1200, bbox_inches='tight')
plt.show()

plt.plot(time,s0, linewidth='3')
plt.xlabel('time (s)') 
plt.ylabel('Concentration (uM)')  
plt.title('o21')
plt.savefig('o212noox.png',format='png',dpi=1200, bbox_inches='tight')
plt.show()
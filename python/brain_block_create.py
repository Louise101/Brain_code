#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 13 15:38:06 2023

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
#fd = open('rhokap_brain_code26.dat', 'rb')
#fd = open('60o21_g21_ox39_thresh560.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
fd.close()
data=datain#np.flip(datain)

#fd = open('600tum_kill_g21_ox39_thresh560.dat', 'rb')
#datain=np.fromfile(file=fd, dtype=np.double).reshape(150,206,236)
#fd.close()
#data_tum=np.flip(datain)

dslice=data[2:154, 115, 88:240]

#block=rho[38:118,90:170,170:]

block_jmean=data[2:154,48:200, 88:240]

jmean_slice=block_jmean[76,:,:]

block_jmean_tran=np.transpose(block_jmean,(1,0,2))

block_slice_jmean=block_jmean_tran[:,76,:]

print(len(block_jmean[0][0]))

block_jmean_flat=block_jmean_tran.flatten()
f=np.savetxt('block_jmean.txt', block_jmean_flat) 


fd = open('rhokap_brain_code26.dat', 'rb')
datain_rho=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
fd.close()
rho=datain_rho#np.flip(datain_rho)

#rho_flat=rho.flatten()
#f=np.savetxt('rho_full.txt', rho_flat)  


fd = open('albedo_brain_code26.dat', 'rb')
datain_alb=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
fd.close()
alb=datain_alb#np.flip(datain_rho)

#alb_flat=alb.flatten()
#f=np.savetxt('alb_full.txt', alb_flat)  

#for i in range (len(rho)):
 # for j in range (len(rho[0])):
  #    for k in range(len(rho[0][0])):
   #      if (rho[i][j][k]==10.001 or rho[i][j][k]==0.020001 or rho[i][j][k]==0.39 or rho[i][j][k]==0. or rho[i][j][k]==0.003):
    #        data[i][j][k]=0.
     #       rho[i][j][k]=0.
      #      alb[i][j][k]=0.
            
#dslice=rho[78,90:170,170:]

#dslice=rho[77,48:200,88:240]

dslice=rho[2:154, 115, 88:240]

#block=rho[38:118,90:170,170:]

block_rho=rho[2:154,48:200, 88:240]

rho_slice=block_rho[76,:,:]

block_rho_tran=np.transpose(block_rho,(1,0,2))

block_slice=block_rho_tran[:,76,:]

print(len(block_rho[0][0]))

block_rho_flat=block_rho_tran.flatten()
f=np.savetxt('block_rho.txt', block_rho_flat)  

#dslice_a=alb[78,90:170,170:]

#block_a=alb[38:118,90:170,170:]



#block_slice_a=block_tran_a[:,40,:]

#print(len(block_tran_a[0][0]))

dslice=alb[2:154, 115, 88:240]

#block=rho[38:118,90:170,170:]

block_alb=alb[2:154,48:200, 88:240]

block_alb_tran=np.transpose(block_alb,(1,0,2))

alb_slice=block_alb[76,:,:]

#block_tran=np.transpose(block,(2,0,1))

#block_slice_a=block_alb[:,40,:]

print(len(block_alb[0][0]))


block_flat_a=block_alb_tran.flatten()
f=np.savetxt('block_alb.txt', block_flat_a)  

plt.imshow(block_slice_jmean)#, norm=color.LogNorm())
plt.xlabel('x (mm)') 
plt.ylabel('y (mm)') 
#plt.title('Jmean')
plt.colorbar()
#plt.savefig('cav_flu.png',format='png',dpi=1200, bbox_inches='tight')
plt.show()
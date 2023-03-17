#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 13 12:08:42 2022

@author: lf58
"""

import matplotlib.pyplot as plt
import numpy as np

fd = open('esc_flur.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double)
fd.close()
data=np.flipud(datain)

fd = open('pl_absor.dat', 'rb')
datain_pl=np.fromfile(file=fd, dtype=np.double)
fd.close()
data_pl=np.flipud(datain_pl)

print(data_pl)
print(data)

tot_dat=sum(data)
print(tot_dat)

prop_dat=[]
for i in range(0,len(data)):
    v=(data[i]/data_pl[i])
    prop_dat.append(v)
    
#print(prop_dat)
z_dep=2. #cm

z= np.arange(0,2.0,2/100)

#print(len(datain),len(z))


e_flu=[]
c3=0.771
k3=0.994
de=0.290

for i in range(len(z)):
    e= c3 * np.exp((-k3 * z[i]/de))
    e_flu.append(e)

    
plt.plot(z,e_flu, label='Predicted')
plt.plot(z,prop_dat, label='MCRT')
plt.legend()
plt.xlabel('depth(cm)')
plt.ylabel('Escape Fluorescence')
plt.show()

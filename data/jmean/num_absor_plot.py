#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 21 10:44:55 2022

@author: lf58
"""

import matplotlib.pyplot as plt
import numpy as np

fd = open('alb_absor.dat', 'rb')
datain_alb=np.fromfile(file=fd, dtype=np.double)
fd.close()
data_alb=np.flipud(datain_alb)

fd = open('pl_absor.dat', 'rb')
datain_pl=np.fromfile(file=fd, dtype=np.double)
fd.close()
data_pl=np.flipud(datain_pl)

fd = open('flu_phot_rel.dat', 'rb')
datain_flu=np.fromfile(file=fd, dtype=np.double)
fd.close()
data_flu=np.flipud(datain_flu)

print(data_pl)

tot_dat=sum(data_flu)
print(tot_dat)   


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
    
    
#plt.plot(z,e_flu)
plt.plot(z,data_flu)
plt.show()

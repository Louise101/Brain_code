#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb  9 16:21:59 2023

@author: lf58
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as color

#fd = open('so_slice24.dat', 'rb')
#fd = open('o23_slice24.dat', 'rb')
#fd = open('o21_slice24.dat', 'rb')
fd = open('tumkill_slice25.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(1200,231,250)
fd.close()
s0=datain#np.flip(datain)



dslice=s0[0,:,:]

plt.imshow(dslice)#, norm=color.LogNorm())
#plt.xlabel('x (mm)') 
#plt.ylabel('y (mm)') 
#plt.title('Jmean')
plt.colorbar()
#plt.savefig('cav_flu.png',format='png',dpi=1200, bbox_inches='tight')
plt.show()
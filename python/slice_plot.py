#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 10 11:39:50 2023

@author: lf58
"""

import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as color

fd = open('so_slice25.dat', 'rb')
#fd = open('o23_slice25.dat', 'rb')
#fd = open('o21_slice25.dat', 'rb')
#fd = open('tumkill_slice25.dat', 'rb')
#fd = open('temp_slice25.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(1200,231,250)
fd.close()
s0=datain#np.flip(datain)

print(len(s0))

dat=s0[0,150,44:68]
depth=np.arange(24)


plt.plot(depth*0.233,dat, linewidth='3', label='O21')
plt.xlabel('Depth') 
plt.ylabel('Concentration')  
#plt.title('o21_35,oxdep')
#plt.legend()
#plt.xticks(fontsize=15)
#plt.yticks(fontsize=15)
#plt.savefig('flu_wall_x_code1_right.png',format='png',dpi=1200, bbox_inches='tight')
plt.show()


print(s0[0,165,44])

data=[]
for j in range(len(s0)):
    data.append(s0[j,165,44])
    
    
time=np.arange(1200)
plt.plot(time,data, linewidth='3', label='O21')
plt.xlabel('Time') 
plt.ylabel('Concentration')  
#plt.title('o21_35,oxdep')
#plt.legend()
#plt.xticks(fontsize=15)
#plt.yticks(fontsize=15)
#plt.savefig('flu_wall_x_code1_right.png',format='png',dpi=1200, bbox_inches='tight')
plt.show()


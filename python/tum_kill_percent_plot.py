#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 13 15:25:53 2022

@author: lf58
"""


import matplotlib.pyplot as plt
import matplotlib.colors as color
import numpy as np
import math as math
import matplotlib.image as mpimg


#fd = open('percent_left_code7.dat', 'rb')
#data_3=np.fromfile(file=fd, dtype=np.double)
#fd.close()

fd = open('percent_left_code17.dat', 'rb')
data_1=np.fromfile(file=fd, dtype=np.double)
fd.close()

fd = open('percent_left_code19.dat', 'rb')
data_3=np.fromfile(file=fd, dtype=np.double)
fd.close()

fd = open('percent_left_code18.dat', 'rb')
data_5=np.fromfile(file=fd, dtype=np.double)
fd.close()

fd = open('percent_left_code16.dat', 'rb')
data_7=np.fromfile(file=fd, dtype=np.double)
fd.close()

fd = open('percent_left_code20.dat', 'rb')
data_10=np.fromfile(file=fd, dtype=np.double)
fd.close()


fd = open('percent_left_code21.dat', 'rb')
data_54w=np.fromfile(file=fd, dtype=np.double)
fd.close()

fd = open('percent_left_code22.dat', 'rb')
data_51w=np.fromfile(file=fd, dtype=np.double)
fd.close()

#fd = open('percent_left_4_7cm_05_100.dat', 'rb')
#data_05=np.fromfile(file=fd, dtype=np.double)
#fd.close()

time=np.arange(1,1201)
time=time/60.

time600=np.arange(1,601)
time600=time600/60.

#time_15=np.arange(1,901)
#time_15=time_15/60.


#plt.plot(time_15,data_2_15, label='2', linewidth=3, color='orange')
#plt.plot(time[:240],data_3[:240], label='Oxygen Depletion', linewidth=3)
#plt.plot(time[:240],data_2[:240], label='Constant Oxygen', linewidth=3)
#plt.plot(time600,data_4, label='Oxygen Depletion 35uM - g=37', linewidth=3)
#plt.plot(time600,data_5, label='no Oxygen Depletion 35uM - g=37', linewidth=3)
#plt.plot(time[:600],data_7[:600], label='Oxygen Depletion, init ox=40um, g=10', linewidth=3)
#plt.plot(time600[:507],data_6[:507], label='Constant Oxygen', linewidth=3)

#plt.plot(time,data_1, label='1',linewidth=3)
#plt.plot(time,data_05, label='0.5',linewidth=3)

#plt.plot(time[:558],data_1[:558], label='1uM', linewidth=3)
#plt.plot(time[:558],data_3[:558], label='3uM', linewidth=3)
plt.plot(time[:558],data_5[:558], label='2W', linewidth=3)
#plt.plot(time[:558],data_7[:558], label='7uM', linewidth=3)
#plt.plot(time[:558],data_10[:558], label='10uM', linewidth=3)

plt.plot(time[:558],data_54w[:558], label='4W', linewidth=3)
plt.plot(time[:558],data_51w[:558], label='1W', linewidth=3)


#plt.yscale('log')
plt.xlabel('Time (Min)') 
plt.ylabel('Percentage of Tumour Remaining')  
plt.title('Percentage of Tumour Remaining over Time')
plt.legend()
plt.xticks(fontsize=18)
plt.yticks(fontsize=18)
plt.savefig('percent_left_power.png',format='png',dpi=1200, bbox_inches='tight')
plt.show()

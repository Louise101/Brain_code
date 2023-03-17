#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 10 11:43:16 2022

@author: lf58
"""

import matplotlib.pyplot as plt
import numpy as np
import math as math
import csv

fd = open('jmean-wang_linear.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(101,100,100)
fd.close()
data=np.flip(datain)



data_x=[]
data_y=[]
data_z=[]

depth_x=[]
depth_y=[]
depth_z=[]

for k in range(len(data)):
    depth_z.append(k*(2/100.))
    data_z.append(data[k][40][50])

for j in range(len(data[0])):
    depth_y.append(j)
    data_y.append(data[50][j][50])
    
for i in range(len(data[0][0])):
    depth_x.append(i*(2/100)-1)
    data_x.append(data[40][50][i])

#for k in range(len(data)):
 #   for j in range(len(data[0])):
  #      for i in range(len(data[0][0])):

   #         data_x.append(data[91][35][i])
    #        data_y.append(data[91][j][58])
     #       data_z.append(data[k][35][58])
     
with open('fluence_wang_linear.csv', newline='')as f:
   reader=csv.reader(f)
   lopez_dat=list(reader)
    
lopez_dep=[]
lopez_flu=[]

for i in range(len(lopez_dat)):
    v=float(lopez_dat[i][0])/10.
    lopez_dep.append(v)
    
for i in range(len(lopez_dat)):
    v=float(lopez_dat[i][1])
    lopez_flu.append(v)           
      
            


    
print(data_z)
    


#0.07cm per voxel
#so if cavity wall at voxel 79 - 10% is 0.25cm and 1% 0.6cm, 33% is at 0.07cm
  
plt.plot(depth_x[50:],data_x[50:],label=('MCRT data'))
#plt.plot(depth_z,data_z,label=('MCRT data'))
plt.plot(lopez_dep, lopez_flu,label=('Wang data'))
plt.legend()
#plt.plot(depth_z[58:69],data_z[58:69])
#plt.yscale('log')
plt.xlabel('x-depth (cm)') 
plt.ylabel('Fluence Rate') 
plt.title('Fluence Rate as a Function of distace from Linear Source')
plt.savefig('linear_source_mouse3_fluence.png',format='png',dpi=1200, bbox_inches='tight')
plt.show()
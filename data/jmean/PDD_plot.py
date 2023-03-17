#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 11 15:51:03 2022

@author: lf58
"""

import matplotlib.pyplot as plt
import numpy as np
import glob
from tkinter import Tcl

c1630=6.27
c2630=1.18
k1630=1.00
k2630=14.4
delta630=0.261
E630=0.0265

irr=0.1
Co=0.2
B630=105.

wl=630*10**(-9)
c=3*10**8
h=6.626*10**-34
b=wl/(h*c)

data630= np.loadtxt('validation-630_1.dat', unpack=True)

x630=(data630[0])
y630=data630[1]


def Yex_630(z):
    return irr*((c1630*np.exp((-k1630*z)/delta630))-(c2630*np.exp((-k2630*z)/delta630)))

def C630(z,t):
    return Co*np.exp((-Yex_630(z)*t)/B630)

def PDD(z,t):
    return b*E630*Co*B630*(1-np.exp(-Yex_630(z)*t/B630))

time630_dat=np.arange(0,153)    
time_tau630_dat=[]
for i in range (len(time630_dat)):
    v=time630_dat[i]/(B630/irr)
    time_tau630_dat.append(v)
    
depth=np.arange(0,2,0.02)


#************* import pdd data***************
file_list=glob.glob('*min_pdd.dat')
#print(file_list)

ordered_files=Tcl().call('lsort','-dict', file_list)
print(ordered_files)

n = 100
data_pdd=[]
for f in ordered_files:
    fd=open(f,'rb')
    datain=np.fromfile(file=fd, dtype=np.double).reshape(n,n,n)
    fd.close()
    flipdata=np.flipud(datain)
    data_pdd.append(flipdata)
    
pdd_1min=data_pdd[0]
pdd_2min=data_pdd[1]
pdd_3min=data_pdd[2]
pdd_4min=data_pdd[3]
pdd_5min=data_pdd[4]
pdd_6min=data_pdd[5]
pdd_7min=data_pdd[6]
pdd_8min=data_pdd[7]
pdd_9min=data_pdd[8]
pdd_10min=data_pdd[9]
pdd_20min=data_pdd[10]
pdd_30min=data_pdd[11]
pdd_40min=data_pdd[12]
pdd_50min=data_pdd[13]
pdd_60min=data_pdd[14]
  
  

#fd = open('128pdd.dat', 'rb')
#datain=np.fromfile(file=fd, dtype=np.double).reshape(n,n,n)
#fd.close()
#data_pdd=np.flipud(datain)



#dslice=data_pdd[:n,0:n,50] #use data for voxel number axis and datain for axis produced in python like depth. 

#for i in range(znumber):
#	idepth=(i*zvoxel_length) #cm
 #   depth.append(idepth)
	

#image=plt.imshow(dslice,extent=[min(depth),max(depth),min(depth),max(depth)])
#plt.imshow(dslice)
#plt.xlabel('x (cm)') 
#plt.ylabel('z (cm)') 
#plt.title('skin layers')
#plt.colorbar()
#plt.show()


#************ 1min *************************************
pdd_lay_1min=[]
for k in range(len(pdd_1min[0][0])):
    v=0.
    for i in range(len(pdd_1min)):
        for j in range(len(pdd_1min[0])):
            v=v+pdd_1min[k][i][j]
            
    pdd_lay_1min.append(v) 

for i in range(len(pdd_lay_1min)):
    pdd_lay_1min[i]=pdd_lay_1min[i]/(n*n)
    

   
#**************2 min *********************************

pdd_lay_2min=[]
for k in range(len(pdd_2min[0][0])):
    v=0.
    for i in range(len(pdd_2min)):
        for j in range(len(pdd_2min[0])):
            v=v+pdd_2min[k][i][j]
            
    pdd_lay_2min.append(v) 

for i in range(len(pdd_lay_2min)):
    pdd_lay_2min[i]=pdd_lay_2min[i]/(n*n)
    
#**************3 min *********************************

pdd_lay_3min=[]
for k in range(len(pdd_3min[0][0])):
    v=0.
    for i in range(len(pdd_3min)):
        for j in range(len(pdd_3min[0])):
            v=v+pdd_3min[k][i][j]
            
    pdd_lay_3min.append(v) 

for i in range(len(pdd_lay_3min)):
    pdd_lay_3min[i]=pdd_lay_3min[i]/(n*n)
    
#**************4 min *********************************

pdd_lay_4min=[]
for k in range(len(pdd_4min[0][0])):
    v=0.
    for i in range(len(pdd_4min)):
        for j in range(len(pdd_4min[0])):
            v=v+pdd_4min[k][i][j]
            
    pdd_lay_4min.append(v) 

for i in range(len(pdd_lay_4min)):
    pdd_lay_4min[i]=pdd_lay_4min[i]/(n*n)
    
#**************5 min *********************************

pdd_lay_5min=[]
for k in range(len(pdd_5min[0][0])):
    v=0.
    for i in range(len(pdd_5min)):
        for j in range(len(pdd_5min[0])):
            v=v+pdd_5min[k][i][j]
            
    pdd_lay_5min.append(v) 

for i in range(len(pdd_lay_5min)):
    pdd_lay_5min[i]=pdd_lay_5min[i]/(n*n)
    
#**************6 min *********************************

pdd_lay_6min=[]
for k in range(len(pdd_6min[0][0])):
    v=0.
    for i in range(len(pdd_6min)):
        for j in range(len(pdd_6min[0])):
            v=v+pdd_6min[k][i][j]
            
    pdd_lay_6min.append(v) 

for i in range(len(pdd_lay_6min)):
    pdd_lay_6min[i]=pdd_lay_6min[i]/(n*n)
    
#**************7 min *********************************

pdd_lay_7min=[]
for k in range(len(pdd_7min[0][0])):
    v=0.
    for i in range(len(pdd_7min)):
        for j in range(len(pdd_7min[0])):
            v=v+pdd_7min[k][i][j]
            
    pdd_lay_7min.append(v) 

for i in range(len(pdd_lay_7min)):
    pdd_lay_7min[i]=pdd_lay_7min[i]/(n*n)
    
#**************8 min *********************************

pdd_lay_8min=[]
for k in range(len(pdd_8min[0][0])):
    v=0.
    for i in range(len(pdd_8min)):
        for j in range(len(pdd_8min[0])):
            v=v+pdd_8min[k][i][j]
            
    pdd_lay_8min.append(v) 

for i in range(len(pdd_lay_8min)):
    pdd_lay_8min[i]=pdd_lay_8min[i]/(n*n)
    
#**************9 min *********************************

pdd_lay_9min=[]
for k in range(len(pdd_9min[0][0])):
    v=0.
    for i in range(len(pdd_9min)):
        for j in range(len(pdd_9min[0])):
            v=v+pdd_9min[k][i][j]
            
    pdd_lay_9min.append(v) 

for i in range(len(pdd_lay_9min)):
    pdd_lay_9min[i]=pdd_lay_9min[i]/(n*n)
    
  
#**************10 min ************************
pdd_lay_10min=[]
for k in range(len(pdd_10min[0][0])):
    v=0.
    for i in range(len(pdd_10min)):
        for j in range(len(pdd_10min[0])):
            v=v+pdd_10min[k][i][j]
            
    pdd_lay_10min.append(v) 

for i in range(len(pdd_lay_10min)):
    pdd_lay_10min[i]=pdd_lay_10min[i]/(n*n)
    
#**************20 min ************************
pdd_lay_20min=[]
for k in range(len(pdd_20min[0][0])):
    v=0.
    for i in range(len(pdd_20min)):
        for j in range(len(pdd_20min[0])):
            v=v+pdd_20min[k][i][j]
            
    pdd_lay_20min.append(v) 

for i in range(len(pdd_lay_20min)):
    pdd_lay_20min[i]=pdd_lay_20min[i]/(n*n)
    
#**************30 min ************************
pdd_lay_30min=[]
for k in range(len(pdd_30min[0][0])):
    v=0.
    for i in range(len(pdd_30min)):
        for j in range(len(pdd_30min[0])):
            v=v+pdd_30min[k][i][j]
            
    pdd_lay_30min.append(v) 

for i in range(len(pdd_lay_30min)):
    pdd_lay_30min[i]=pdd_lay_30min[i]/(n*n)
    
#**************40 min ************************
pdd_lay_40min=[]
for k in range(len(pdd_40min[0][0])):
    v=0.
    for i in range(len(pdd_40min)):
        for j in range(len(pdd_40min[0])):
            v=v+pdd_40min[k][i][j]
            
    pdd_lay_40min.append(v) 

for i in range(len(pdd_lay_40min)):
    pdd_lay_40min[i]=pdd_lay_40min[i]/(n*n)
    
#**************50 min ************************
pdd_lay_50min=[]
for k in range(len(pdd_50min[0][0])):
    v=0.
    for i in range(len(pdd_50min)):
        for j in range(len(pdd_50min[0])):
            v=v+pdd_50min[k][i][j]
            
    pdd_lay_50min.append(v) 

for i in range(len(pdd_lay_50min)):
    pdd_lay_50min[i]=pdd_lay_50min[i]/(n*n)
    
#**************60 min ************************
pdd_lay_60min=[]
for k in range(len(pdd_60min[0][0])):
    v=0.
    for i in range(len(pdd_60min)):
        for j in range(len(pdd_60min[0])):
            v=v+pdd_60min[k][i][j]
            
    pdd_lay_60min.append(v) 

for i in range(len(pdd_lay_60min)):
    pdd_lay_60min[i]=pdd_lay_60min[i]/(n*n)
   
#*************************Plot create*********************************************    
plt.plot(depth,PDD(depth,60.))#, label='1 min')
plt.plot(depth,PDD(depth,120.), label='Predicted Model')
plt.plot(depth,PDD(depth,180.))#, label='3 min predicted')
plt.plot(depth,PDD(depth,240.))#, label='4 min predicted')
plt.plot(depth,PDD(depth,300.))#, label='5 min predicted')
plt.plot(depth,PDD(depth,360.))#, label='6 min predicted')
plt.plot(depth,PDD(depth,420.))#, label='7 min predicted')
plt.plot(depth,PDD(depth,480.))#, label='8 min predicted')
plt.plot(depth,PDD(depth,540.))#, label='9 min predicted')
plt.plot(depth,PDD(depth,600.))#, label='10 min')
plt.plot(depth,PDD(depth,1200.))#, label='20 min predicted')
plt.plot(depth,PDD(depth,1800.))#, label='30 min predicted')
plt.plot(depth,PDD(depth,2400.))#, label='40 min predicted')
plt.plot(depth,PDD(depth,3000.))#, label='50 min predicted')
plt.plot(depth,PDD(depth,3600.))#, label='60 min')

plt.plot(depth,pdd_lay_1min,'k--', label='MCRT')
plt.plot(depth,pdd_lay_2min,'k--')#, label='2 min MCRT')
plt.plot(depth,pdd_lay_3min,'k--')#, label='3 min MCRT')
plt.plot(depth,pdd_lay_4min,'k--')#, label='4 min MCRT')
plt.plot(depth,pdd_lay_5min,'k--')#, label='5 min MCRT')
plt.plot(depth,pdd_lay_6min,'k--')#, label='6 min MCRT')
plt.plot(depth,pdd_lay_7min,'k--')#, label='7 min MCRT')
plt.plot(depth,pdd_lay_8min,'k--')#, label='8 min MCRT')
plt.plot(depth,pdd_lay_9min,'k--')#, label='9 min MCRT')
plt.plot(depth,pdd_lay_10min,'k--')#, label='10 min MCRT')
plt.plot(depth,pdd_lay_20min,'k--')#, label='20 min MCRT')
plt.plot(depth,pdd_lay_30min,'k--')#, label='30 min MCRT')
plt.plot(depth,pdd_lay_40min,'k--')#, label='40 min MCRT')
plt.plot(depth,pdd_lay_50min,'k--')#, label='50 min MCRT')
plt.plot(depth,pdd_lay_60min,'k--')#, label='60 min MCRT')



#plt.yscale('log')
plt.legend()
plt.xlabel('Depth(cm)')
plt.ylabel('PDD (photons/cm^3)')
plt.show()

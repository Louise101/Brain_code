#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 27 14:12:35 2022

@author: lf58
"""

import matplotlib.pyplot as plt
import numpy as np
import glob
from tkinter import Tcl

#code produces theoratical ovserved fluorescence based on eq 4 in Jaques 


irr=0.1
c1420=5.76
c2420=1.31
k1420=1.00
k2420=10.2
delta420=0.047
E420=0.105


c1630=6.27
c2630=1.18
k1630=1.00
k2630=14.4
delta630=0.261
E630=0.0265

Co=0.2
B630=105.#14.
ua630=0.06#0.105*0.2
ua420=1.1 # from fig 1.4 in Campbell thesis
B420=105.#5.*(B630/3.6)#(B630*ua630*630.)/(420.*ua420) #from chap 3 in campbell thesis
print(B420)
fy=1. #fluorescence yield
f=1. #collection efficiency

data420= np.loadtxt('validation-420_1.dat', unpack=True)
data630= np.loadtxt('validation-630_1.dat', unpack=True)

x420=data420[0]
y420=data420[1]
x630=(data630[0])
y630=data630[1]


#theoretical excitiation fluence rate Yex(z)

def Yex_420(z):
    return irr*((c1420*np.exp((-k1420*z)/delta420))-(c2420*np.exp((-k2420*z)/delta420)))

def Yex_630(z):
    return irr*((c1630*np.exp((-k1630*z)/delta630))-(c2630*np.exp((-k2630*z)/delta630)))

#PPIX concentration

def C420(z,t):
    return Co*np.exp((-Yex_630(z)*t)/B420) #needs fluence of the photobleaching light which is 630nm (sec 3.2 on Jaques paper)

def C630(z,t):
    return Co*np.exp((-Yex_630(z)*t)/B630)

time420=np.arange(0,2100) #running for 35 minutes
time_tau420=[]
for i in range (len(time420)):
    v=time420[i]/(B420/irr)
    time_tau420.append(v)

#time630=np.arange(0,153)   
time630=[]
t=0
for i in range(3600):
    if i < 3601:
        t=t+1
        time630.append(t)
    else:
        t=t+60
        time630.append(t)


 
time_tau630=[]

for i in range (len(time630)):
    v=time630[i]/(B630/irr)
    time_tau630.append(v)
      
#Escape fluorescence E(z)
c3=0.771
k3=0.994
de=0.290

def E(z):
    return c3 * np.exp((-k3 * z/de))

f630_0=0.    
for j in range(len(x630)):
        f630_0= f630_0 + Yex_630(x630[j])*E630*C630(x630[j],0.)*fy*E(x630[j])*f 
    
F_630=[] 
for i in range(len(time630)):
    f630v=0.
    for j in range(len(x630)):
        f630v= f630v + Yex_630(x630[j])*E630*C630(x630[j],time630[i])*fy*E(x630[j])*f#/(Yex_420*E420*C420[0]*fy*Ez420*f ) 
    F_630.append(f630v/f630_0)
    
f420_0=0.    
for j in range(len(x420)):
        f420_0= f420_0 + Yex_420(x420[j])*E420*C420(x420[j],0.)*fy*E(x420[j])*f 
 
F_420=[] 
for i in range(len(time630)):
    f420v=0.
    for j in range(len(x420)):
        f420v= f420v + (Yex_420(x420[j])*E420*C420(x420[j],time630[i])*fy*E(x420[j])*f)#/(Yex_420(x420[j])*E420*C420(x420[j],0.)*fy*E(x420[j])*f) 
    F_420.append(f420v/f420_0)
    
#******************************* MCRT results plot*********************************

#file_list=glob.glob('*obs_flu.dat')
#ordered_files=Tcl().call('lsort','-dict', file_list)
#print(ordered_files)

fd = open('obs_flur.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double)
fd.close()
data=np.flipud(datain)

fd = open('obs_flur420.dat', 'rb')
datain420=np.fromfile(file=fd, dtype=np.double)
fd.close()
data420=np.flipud(datain420)



norm_dat=[]
for i in range(len(data)):
    v=datain[i]/datain[0]
    norm_dat.append(v)

norm_dat420=[]
for i in range(len(data420)):
    v=datain420[i]/datain420[0]
    norm_dat420.append(v)
    
    
time630_dat=np.arange(0,3600)    
time_tau630_dat=[]
for i in range (len(time630_dat)):
    v=time630_dat[i]/(B630/irr)
    time_tau630_dat.append(v)
  





  
#*************************Plot create*********************************************    
plt.plot(time_tau630,F_630, label='Predicted')
#plt.plot(time_tau630,F_420, label='420nm predicted')
plt.plot(time_tau630,norm_dat,'r--', label='MCRT')
#plt.plot(time_tau630,norm_dat420,'g--', label='MCRT 420nm')
#plt.yscale('log')
plt.legend()
plt.xlabel('Light Exposure (t/tau)')
plt.ylabel('Observed Fluorescence')
plt.show()
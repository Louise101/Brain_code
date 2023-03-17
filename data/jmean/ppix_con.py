#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  8 15:04:52 2022

@author: lf58
"""

import matplotlib.pyplot as plt
import numpy as np


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
B420=5.*(B630/3.6)#(B630*ua630*630.)/(420.*ua420) #from chap 3 in campbell thesis

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

time420=np.arange(0,300) #running for 35 minutes
time_tau420=[]
for i in range (len(time420)):
    v=time420[i]/(B420/irr)
    time_tau420.append(v)

#time630=np.arange(0,153)    

time630=[]
ti=0
for i in range(153):
    if i < 121:
        ti=ti+1
        time630.append(ti)
    else:
        ti=ti+60
        time630.append(ti)
        
time_tau630=[]
for i in range (len(time630)):
    v=time630[i]/(B630/irr)
    time_tau630.append(v)
    
fd = open('con_test.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double)
fd.close()
data=np.flipud(datain)

print(datain)
    
plt.plot(time630,C630(x630[19],time630), label='630nm predicted')
#plt.plot(time_tau420,F_420, label='420nm predicted')
plt.plot(time630,datain, label='MCRT 630nm')
#plt.yscale('log')
plt.legend()
plt.xlabel('Light Exposure (t/tau)')
plt.ylabel('ppix concentration')
plt.show()
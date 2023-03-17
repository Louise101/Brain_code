#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 19 16:56:55 2023

@author: lf58
"""

import matplotlib.pyplot as plt
import numpy as np
#calculations from tissue oxygen supply rate paper zhu et al. 

vz=200.
Rc=2.5
Rt=30. #max - likely to be different as vessels in brain are not evenly spaced but that is whta krough model assumes
lz=200.
q0=26.3

g100=(1200*vz*Rc*(Rc + ((100**2 + q0**2)/(50**2 -q0**2))))/(lz*(Rt+4.2)**2)# - 26.25

g50=(1200*vz*Rc*(Rc + ((50**2 + q0**2)/(50**2 +q0**2))))/(lz*(Rt-4.2)**2)# - 26.25
#g100_0=(1200*vz*Rc*(Rc + ((100**2 + 1.9**2)/(50**2 -1.9**2))))/(lz*(Rt+4.2)**2)# - 26.25



#print(g100,g50)
t=1

td=(t-750)/632.1
g=g100*((0.99*td**4+1.09*td**3+0.05*td**2+0.18*td+0.32)/(td**4+1.16*td**3+0.18*td**2+0.24*td +0.31))

print(g100)


o23_0=83
o23=83
epsil=3.7e-3
jme=300.
s0=7.
beta=11.9
sig=9E-5
delta=33.

diff= (g)*(1-(o23/o23_0)) - (epsil*jme*s0/(o23+beta))*o23#-(maxmet*(o23/(k50+o23)))
#diff= (g50)*(1-(o23/o23_0))-(epsil*jme*s0/(o23+beta))*o23#-(maxmet*(o23/(k50+o23)))

sub=(epsil*jme*s0/(o23+beta))*o23# -(maxmet*(o23/(k50+o23)))

val=g*(1-(o23/o23_0))

print(diff, val, sub)

o21tot=0.
o21=[]
o23_loop=83.
o23_ar=[]
so_ar=[]
#singlet oxygen
for i in range(600):
    
    s0=s0 + (-epsil*sig*jme*(s0+delta)*o23_loop/(o23_loop+beta))*s0
    so_ar.append(s0)
    g=g100*((0.99*i**4+1.09*i**3+0.05*i**2+0.18*i+0.32)/(i**4+1.16*i**3+0.18*i**2+0.24*i +0.31))
    diffloop= g*(1-(o23_loop/o23_0)) - (epsil*jme*s0/(o23_loop+beta))*o23_loop
    o23_loop=80#o23_loop + diffloop
    o23_ar.append(o23_loop)
    
    o21add=(epsil*jme*s0*o23_loop/(o23_loop+beta))
    o21tot = o21tot + o21add
    o21.append(o21tot)
   
o23_loop=40.
o21adda=(epsil*jme*s0*o23_loop/(o23_loop+beta))
print(o21adda)

time=np.arange(600)
#print(so_ar,o23_ar,o21)
#print(o21)
   
#plt.plot(time,so_ar, linewidth='3', label='s0')
#plt.plot(time,o23_ar, linewidth='3', label='O23')
plt.plot(time,o21, linewidth='3', label='O21')
plt.xlabel('Time') 
plt.ylabel('Concentration')  
plt.title('o21_35,oxdep')
#plt.legend()
#plt.xticks(fontsize=15)
#plt.yticks(fontsize=15)
#plt.savefig('flu_wall_x_code1_right.png',format='png',dpi=1200, bbox_inches='tight')
plt.show()
    

s0=35
o23=22
jme=300.
epsil=3.7e-3
o211=(epsil*jme*s0*o23/(o23+beta))

print(o211)
   
   
   
   
   
   
   
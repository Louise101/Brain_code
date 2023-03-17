#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 14:32:42 2023

@author: lf58
"""



import matplotlib.pyplot as plt
import numpy as np
import matplotlib.colors as color

fd = open('rhokap_brain_hires.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
fd.close()
rhokap=np.flip(datain)

rhoslice=rhokap[77,:,:]

fd = open('jmean_brain_g21_ox39_hires_ng.dat', 'rb')
#fd = open('60o21_g21_ox39_thresh560_con023_hires.dat', 'rb')
datain=np.fromfile(file=fd, dtype=np.double).reshape(156,231,250)
fd.close()
jmean=datain#np.flip(datain)

print(np.max(jmean))
print(np.min(jmean))


for i in range (len(rhokap)):
  for j in range (len(rhokap[0])):
      for k in range(len(rhokap[0][0])):
          if (rhokap[i][j][k]==0.38999999999999996 or rhokap[i][j][k]==10.001 or rhokap[i][j][k]==1E-4 + 0.3 or rhokap[i][j][k]==0.or rhokap[i][j][k]== 3.00000000e-03 ):
            jmean[i][j][k]=0.
            
            
 
depth_x=[]
jmean_x=[]
dose_x=[]
             
depth_z=[]
jmean_z=[]
dose_z=[]

depth_y=[]
jmean_y=[]
dose_y=[]

for i in range(len(rhokap[0][0])):
    depth_x.append(i)
    jmean_x.append(rhokap[90][40][i])
    dose_x.append(rhokap[90][40][i]*558)

for k in range(len(rhokap)):
    depth_z.append(k)
    jmean_z.append(rhokap[k][40][56])
    dose_z.append(rhokap[k][40][56]*558)

for j in range(len(rhokap[0])):
    depth_y.append(j)
    jmean_y.append(rhokap[90][j][56])
    dose_y.append(rhokap[90][j][56]*558)
               
#*************** image xy *********************************************
dslice_xy=jmean[77,:,:]
    
centre_y=125#(206-40)*0.7

line_xy=[]
for i in range(len(depth_x)):
    line_xy.append(centre_y) # use y centre
    
image=plt.imshow(dslice_xy,extent=[min(depth_x),max(depth_x),min(depth_y),max(depth_y)],cmap='Reds',norm=color.LogNorm())	


#plt.axis('off')
#image.axes.get_xaxis().set_visible(False)
#image.axes.get_yaxis().set_visible(False)
plt.savefig('jmean_line_xy.png',format='png',dpi=1200, bbox_inches='tight', pad_inches=0)
plt.show('off')

plt.plot(depth_x,line_xy, color='r')
im = plt.imread('jmean_line_xy.png');
plt.xlabel('x (cm)') 
plt.ylabel('y (cm)') 
#plt.title('1 Sec - Heterogeneous')
#plt.colorbar()
plt.imshow(im, interpolation='nearest', zorder=0,extent=[min(depth_x),max(depth_x),min(depth_y),max(depth_y)])
plt.show()

#************** image zx ************************************************
dslice_zx=jmean[:,115,:]

centre_x=125#(150-89)*0.7


line_zx=[]
for i in range(len(depth_z)):
    line_zx.append(centre_x)
    
image=plt.imshow(dslice_zx,extent=[min(depth_x),max(depth_x),min(depth_z),max(depth_z)],norm=color.LogNorm())	
plt.axis('off')
image.axes.get_xaxis().set_visible(False)
image.axes.get_yaxis().set_visible(False)
plt.savefig('jmean_line_zx.png',format='png',dpi=1200, bbox_inches='tight', pad_inches=0)
plt.show('off')

plt.plot(line_zx,depth_z, color='r')
im = plt.imread('jmean_line_zx.png');
plt.xlabel('z (cm)') 
plt.ylabel('x (cm)') 
#plt.title('1 Sec - Heterogeneous')
#plt.colorbar()
plt.imshow(im, interpolation='nearest', zorder=0,extent=[min(depth_x),max(depth_x),min(depth_z),max(depth_z)])
plt.show()

#************** image yz ***********************************************
dslice_y=jmean[:,:,125]
centre_z=77#(236-56)*0.7 
line_zy=[]
for i in range(len(depth_y)):
    line_zy.append(centre_z)

image=plt.imshow(dslice_y,extent=[min(depth_y),max(depth_y),min(depth_z),max(depth_z)],norm=color.LogNorm())	
plt.axis('off')
image.axes.get_xaxis().set_visible(False)
image.axes.get_yaxis().set_visible(False)
plt.savefig('jmean_line_y.png',format='png',dpi=1200, bbox_inches='tight', pad_inches=0)
plt.show('off')

plt.plot(depth_y,line_zy, color='r')
im = plt.imread('jmean_line_y.png');
plt.xlabel('z (cm)') 
plt.ylabel('x (cm)') 
#plt.title('1 Sec - Heterogeneous')
#plt.colorbar()
plt.imshow(im, interpolation='nearest', zorder=0,extent=[min(depth_y),max(depth_y),min(depth_z),max(depth_z)])
plt.show()    

    




#*************** Jmean / Fluence rate *******************************************************
#fd = open('jmean.dat', 'rb')
#datain=np.fromfile(file=fd, dtype=np.double).reshape(150,206,236)
#fd.close()
#jmean=np.flip(datain)


#make non brain/tumour equal 0. 
#for i in range (len(rhokap)):
 #   for j in range (len(rhokap[0])):
  #      for k in range(len(rhokap[0][0])):
   #         if (rhokap[i][j][k]==0.38999999999999996 or rhokap[i][j][k]==10.001):
    #           jmean[i][j][k]=0.


depth_x=[]
jmean_x=[]
dose_x=[]
rho_x=[]
             
depth_z=[]
jmean_z=[]
dose_z=[]

depth_y=[]
jmean_y=[]
dose_y=[]

for k in range(len(jmean)):
    depth_z.append(k*0.233)
    jmean_z.append(jmean[k][115][125])
    dose_z.append(jmean[k][115][125])#*(558/1000)) #J/cm2

for j in range(len(jmean[0])):
    depth_y.append(j*0.233)
    jmean_y.append(jmean[77][j][125])
    dose_y.append(jmean[77][j][125] )#*(558/1000))
    
for i in range(len(jmean[0][0])):
    depth_x.append(i*0.233)
    jmean_x.append((jmean[77][115][i]))
    dose_x.append((jmean[77][115][i]))#*(558/1000))
    rho_x.append((rhokap[77][115][i]))
    
    
print('rhokap',rho_x[85])
    
    
jmean_x_depth=[]
jmean_x_flu=[]
jmean_x_wall=[]
jmean_x_25=[]
wall=58.8 - 58.8
flu=60.9 - 58.8
for i in range(85,105):
    jmean_x_depth.append(25) # lille 25 threshold for PDT success
    jmean_x_flu.append(1.5)
    jmean_x_25.append(7.1)
    jmean_x_wall.append(wall)
    depth_x[i]=depth_x[i]-58.8
    
#for i in range(len(depth_x)):
 #   depth_x[i]=depth_x[i]-58.8

#89,40,56
#rint(depth_x[95:106])

#plt.plot(depth_x[94:105],jmean_x[94:105])
#plt.plot(depth_x[86:95],jmean_x[86:95])
#plt.plot(depth_x[84:98],jmean_x[84:98])
plt.plot(depth_x,dose_x, linewidth='3', label='Heterogeneous')

#plt.plot(depth_x[85:105],jmean_x_depth, linestyle='dashed', color='black')
#plt.plot(jmean_x_flu,jmean_x[85:105])
#plt.plot(jmean_x_25,jmean_x[85:105])
#plt.plot(jmean_x_wall,jmean_x[84:105])
#plt.yscale('log')
#plt.xlabel('x-depth (mm)') 
#plt.ylabel('Fluence Rate')  
#plt.title('x depth fluence rate ')
#plt.legend()
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.savefig('flu_wall_x_4cm_3_2mW_100_het.png',format='png',dpi=1200, bbox_inches='tight')
plt.show()

x_depth_flu=60.9-58.8
x_wall_flu=jmean_x[84]                
print('x_depth=',x_depth_flu, 'wall flu=', x_wall_flu)

 

jmean_y_depth=[]
jmean_y_flu=[]
jmean_y_wall=[]
jmean_y_25=[]
wall=54.6 - 54.6
flu=55.75 - 54.6
for i in range(78,99):
    jmean_y_depth.append(25)
    jmean_y_flu.append(0.7)
    jmean_y_wall.append(wall)
    jmean_y_25.append(5)
    depth_y[i]=depth_y[i]-54.6

#'plt.plot(depth_y[75:83],jmean_y[75:83])
#plt.plot(depth_y[78:87],jmean_y[78:87])
plt.plot(depth_y[78:99],dose_y[78:99],linewidth='3')

plt.plot(depth_y[78:99],jmean_y_depth, linestyle='dashed', color='black')
plt.plot(jmean_y_25,jmean_y[78:99])
plt.plot(jmean_y_flu,jmean_y[78:99])
#plt.plot(jmean_y_wall,jmean_y[78:100])
#plt.yscale('log')
#plt.xlabel('y-depth (mm)') 
#plt.ylabel('Fluence Rate')  
#plt.title('y depth fluence rate')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.savefig('flu_wall_y_4cm_3_2mW_100_het.png',format='png',dpi=1200, bbox_inches='tight')
plt.show()

y_depth_flu=55.75-54.6
y_wall_flu=jmean_y[78]
print('y_depth=',y_depth_flu, 'wall flu=', y_wall_flu)

jmean_z_depth=[]
jmean_z_flu=[]
jmean_z_wall=[]
jmean_z_25=[]
wall=77-77
flu=79.5-77
for i in range(111,131):
    jmean_z_depth.append(25)
    jmean_z_flu.append(1.8)
    jmean_z_wall.append(wall)
    jmean_z_25.append(5.8)
    depth_z[i]=depth_z[i]-77

#'plt.plot(depth_y[75:83],jmean_y[75:83])
#plt.plot(depth_z[110:125],jmean_z[110:125])
plt.plot(depth_z[111:131],dose_z[111:131],linewidth='3')

plt.plot(depth_z[111:131],jmean_z_depth,linestyle='dashed', color='black')
plt.plot(jmean_z_25,jmean_z[111:131])
plt.plot(jmean_z_flu,jmean_z[111:131])
#plt.plot(jmean_z_wall,jmean_z[110:125])
#plt.yscale('log')
#plt.xlabel('z-depth (mm)') 
#plt.ylabel('Fluence Rate')  
#plt.title('z depth fluence rate')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)
plt.savefig('flu_wall_z_4cm_3_2mW_100_het.png',format='png',dpi=1200, bbox_inches='tight')
plt.show()

z_depth_flu=79.5-77
z_wall_flu=jmean_z[110]
print('z_depth=',z_depth_flu, 'wall flu=', z_wall_flu)
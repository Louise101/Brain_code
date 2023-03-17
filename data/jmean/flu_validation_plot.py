import matplotlib.pyplot as plt
import numpy as np
import math as math
import matplotlib.image as mpimg
import glob

depth=[]
avhole420=[]
avhole630=[]
data=[]
holes420=[]
holes630=[]
holelay420=[]
holelay630=[]
div420=[]
div630=[]

n = 200
#h=6.62607004*10**(-30) #cm^2 kg/s
#c=2.99792458*10**(10) #cm/s
#nphotons=2000000 
#ri=1.38
#P=1

c1420=5.76
c2420=1.31
k1420=1.00
k2420=10.2
delta420=0.047


c1630=6.27
c2630=1.18
k1630=1.00
k2630=14.4
delta630=0.261


zlength= 2. #cm
znumber=100
zvoxel_length=zlength/(znumber) #cm
vox_vol=(zvoxel_length)**3

#L=P*(zlength**2)

#incIr420=((nphotons*h*c)/420*10**(-7))/(zlength**2)
#incIr630=((nphotons*h*c)/630*10**(-7))/(zlength**2)

file_list=glob.glob('validation*.dat')
print(file_list)



#for f in file_list:
#	fd = open('validation-420.dat', 'rb')
#	datain=np.fromfile(file=fd, dtype=np.double)#.reshape(n,n,n)
#	fd.close()
	#flipdata=np.flipud(datain)
	#data.append(datain)

#x420, y420 = np.loadtxt('validation-420.dat', unpack=True)
#x630, y630 = np.loadtxt('validation-630.dat', unpack=True)

#fd = open('validation-420.dat', 'rb')
#datain=np.fromfile(file=fd, dtype=np.double)#.reshape(n,n,n)
#fd.close()
#flipdata=np.flipud(datain)
#data.append(datain)

data420= np.loadtxt('validation-420_1.dat', unpack=True)
data630= np.loadtxt('validation-630_1.dat', unpack=True)


#print(data420[1])

x420=data420[0]
y420=data420[1]
x630=(data630[0])
y630=data630[1]

#print(y420)
#data420=data[0]
#data630=data[1]
    

# average all holes
#for i in range(n):
#	for j in range(n):
#		holev=data420[0:n,i,j]
#		holes420.append(holev)
        
	

#for i in range(n):
#	for j in range(n*n):
#		avholev=holes420[j][i]
#		holelay420.append(avholev)
		
#holelay420a=np.reshape(holelay420,(n,n*n))


#for i in range(n):
#	val=np.sum(holelay420a[i])/(n*n)
#	avhole420.append(val)
	
#for i in range(n):
#	for j in range(n):
#		holev=data630[0:n,i,j]
#		holes630.append(holev)
        

#for i in range(n):
#	for j in range(n*n):
#		avholev=holes630[j][i]
#		holelay630.append(avholev)
		
#holelay630a=np.reshape(holelay630,(n,n*n))

#for i in range(n):
#	val=np.sum(holelay630a[i])/(n*n)
#	avhole630.append(val)
	
#Do maths to calculate fluence rate from path data
#avhole420s=np.asarray(avhole420)
#avhole630s=np.asarray(avhole630)

#hole420s=np.asarray(holes420)
#hole630s=np.asarray(holes630)

#j420=(avhole420s*L/(nphotons*vox_vol))
#j630=(avhole630s*L/(nphotons*vox_vol))

#j420=avhole420s
#j630=avhole630s

#print(L/(nphotons*vox_vol))
#use data for voxel number axis and datain for axis produced in python like depth. 	
#for i in range(znumber):
#	idepth=(i*zvoxel_length) #cm
#	depth.append(idepth)

#depths=np.asarray(depth)

#plot functions
fl420=(c1420*np.exp((-k1420*x420)/delta420)-c2420*np.exp((-k2420*x420)/delta420))

fl630=(c1630*np.exp((-k1630*x630)/delta630)-c2630*np.exp((-k2630*x630)/delta630))



#for i in range(n*n):
	#plt.plot(depths, holes420[i])
	
#for i in range(n*n):
	#plt.plot(depths, holes630[i])

#for i in range (len(depths)):
 #   div420v=fl420[i]/j420[i]
  #  div420.append(div420v)
    
#for i in range (len(depths)):
 #   div630v=fl630[i]/j630[i]
  #  div630.append(div630v)
    
#print(div630)
#print (j630)
    
#plt.plot(depths, div420, label=('420nm'))
#plt.plot(depths, div630, label=('630nm'))
	
#plt.plot(depths, avhole420s, label=('420nm average'))
#plt.plot(depths, avhole630s, label=('630nm average'))

	
plt.plot(x420, fl420, label=('420nm predicted'))

plt.plot(x630, fl630, label=('630nm predicted'))
	
plt.plot(x420, y420, label=('420nm MCRT'), linestyle='dashed')

plt.plot(x630, y630/0.1, label=('630nm MCRT'),linestyle='dashed') #normalised by incident irradience 0.1
	
plt.legend()
plt.title('Jacques Fluence Validation')
plt.xlabel('depth(cm)')
plt.ylabel('Normalised Fluence Rate')
plt.show()

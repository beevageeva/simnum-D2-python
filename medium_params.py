import numpy as np
import math
from constants import gamma, z0, zf
from common import getArrayZShape
from sys import exit


p00 = 10.0
v00 = 0.0

#mediumType = "homog"
mediumType = "inhomog"  #variable density rho00 to test with wave packet

if(mediumType=="homog"):
	#rho00 = 2.0
	#rho00 = 0.5
	rho00 = 1.0
	cs00 = math.sqrt(gamma * p00 / rho00)

elif(mediumType=="inhomog"):

	#rhoType = 1
	rhoType = 4
	#the following TODO not used	 (grad cs constante)
	#rhoType = 3
	#rhoType = 2


	#define densargfunc here!
	#like argFunc:
#	densargFunc = argFunc
	#diagonal
#	k1 = 1.0
#	k2 = 1.0
#	modk = (k1**2 + k2**2)**0.5
#	densargFunc = lambda x: k1/modk * x[0] + k2/modk * x[1]	
	#vertical
	densargFunc = lambda z: z[0] 

	if rhoType == 1:
		inhomogSubtype =1
		#inhomogSubtype =2
		if inhomogSubtype == 1:
			rho00 = 1.0
			rho01 = 0.01
		else:
			rho00 = 0.01
			rho01 = 1.0

		#rho00 = 0.3  #second exp of inhom
		#rho01 = 1.2#second exp of inhom
		ze = [0.5*(z0[0] + zf[0]),0.5*(z0[1] + zf[1])]
		#ze = [z0[0] + 0.7*(zf[0] - z0[0]),z0[1] + 0.7*(zf[1] - z0[1])]#second exp of inhom
		we = 5.0
		#we = 0.5#second exp of inhom
		#I have to apply func (argument function) before applying tanh: see initcond_soundwave
		densFunc = lambda z: 1 + np.tanh(z /we)
	
		
	
		
		def rho0(z):
			from common import getArrayZShape
			if hasattr(z, "__len__") and hasattr(z[0], "__len__")  :	
				zarg = z - getArrayZShape(ze[0], ze[1], len(z[0]))
			else:
				zarg = [z[0] - ze[0], z[1] - ze[1]]
			return rho00 + 0.5 * (rho01-rho00) * densFunc(densargFunc(zarg))

		#TODO why was this here?
		#cs00 = lambda z: np.sqrt(gamma * p00 / densFunc(densargFunc(z)))


	elif rhoType == 4:
		ze = [0.02*(z0[0] + zf[0]),0.02*(z0[1] + zf[1])]
		rho00 = 0.05
		rho01 = 0.8
		#inversion termica
		p00 = 1.01325e5
		T01 = 283.0  #10
		T00 = 323.0  #50
		#mm = 28.97
		#R = 8.314472
		#rho00 = (p00 * mm) / (T00 * R)
		#rho01 = (p00 * mm) / (T01 * R)
		
		Rsp = 287.058
		rho00 = p00/ (T00 * Rsp)
		rho01 = p00/ (T01 * Rsp)
		#end inversion termica
		we = 0.07
		#we = 0.5#second exp of inhom
		#I have to apply func (argument function) before applying tanh: see initcond_soundwave
		densFunc = lambda z: 1 + np.tanh(z /we)
		def rho0(z):
			if hasattr(z, "__len__") and hasattr(z[0], "__len__")  :	
				n1 = z[0].shape[0]
				n2 = z[0].shape[1]
				#print("n1=%d, n2=%d" % (n1,n2))
				mid = int(n1/2)
				res = np.zeros((n1,n2))
				for i in range(mid): 
					for j in range(n2):
						res[i,j] =  rho00 + 0.5 * (rho01-rho00) * densFunc(densargFunc([z[0][i,j] - ze[0], z[1][i,j] - ze[1]]))
				for i in range(mid): 
					for j in range(n2):
						#res[i,j] =  rho01 + 0.5 * (rho00-rho01) * densFunc(densargFunc([z[0][i,j] -0.5*(zf[0]-z0[0]) - 1.0 + ze[0], z[1][i,j] -0.5 *(zf[1]-z0[1]) - 1.0 + ze[1] ]))
						#put reverse
						res[mid+i,j] = res[mid-i-1,j]
				#print("temperature at 0 %e" % ((p00 * mm)/(res[0,0] * R)))	
				#print("temperature at middle %e" % ((p00 * mm)/(res[mid,0] * R)))

				#print always	
				#print("temperature at 0 %e" % (p00 /(res[0,0] * Rsp)))	
				#print("temperature at middle %e" % (p00 /(res[mid,0] * Rsp)))	
				return res
			else:
				if( densargFunc([z[0] - z0[0], z[1] - z0[1]]) < 0.5 *  densargFunc([zf[0] - z0[0], zf[1] - z0[1]] )):
					return rho00 + 0.5 * (rho01-rho00) * densFunc(densargFunc([z[0] - ze[0], z[1] - ze[1]]))
				return rho01 + 0.5 * (rho00-rho01) * densFunc(densargFunc([z[0] - ze[0], z[1] - ze[1]]))

	
	#TODO the following not used
	elif rhoType == 2:
		#c of form -xi + A
		rho00 = 1.0
		a = 3.0
		A = np.sqrt(gamma * p00 / rho00)
		def rho0(z):
			return  p00*gamma/(A + a * densargFunc(z))**2 
		cs00 =  lambda z: A + a * densargFunc(z)

	
	elif rhoType == 3:
		#c of form -xi + A
		#inversion termica
		p00 = 1.01325e5
		#mm = 1.660538921e-27 this is atomic weight!
		mm = 28.97
		T00 = 293.0  #20
		R = 8.314472
		rho00 = (p00 * mm) / (T00 * R)
		#end inversion termica
		print("rho00 in rhoType = 3 = %e" % rho00)
		a = 0.4
		A = np.sqrt(gamma * p00 / rho00)
		from constants import z0,zf
		def rho0(z):
			zarg = z[0]
			if(not hasattr(zarg,"__len__")):
				if(zarg-z0[0]<0.5*(zf[0]-z0[0])):
					return p00*gamma/(A - a * zarg)**2  
				return  p00*gamma/(A + a * zarg)**2
			n1 = zarg.shape[0]
			n2 = zarg.shape[1]
			middle = int(n1/2)
			res = np.zeros((n1,n2))
			for i in range(middle):
				for j in range(n2):
					res[i,j] = p00*gamma/(A - a * zarg[i,j])**2

			for i in range(middle): 
				for j in range(n2):
						#res[i,j] =  rho01 + 0.5 * (rho00-rho01) * densFunc(densargFunc([z[0][i,j] -0.5*(zf[0]-z0[0]) - 1.0 + ze[0], z[1][i,j] -0.5 *(zf[1]-z0[1]) - 1.0 + ze[1] ]))
						#put reverse
					res[middle+i,j] = res[middle-i-1,j]

			print("temperature at 0 %e" % ((p00 * mm)/(res[0,0] * R)))	
			print("temperature at middle %e" % ((p00 * mm)/(res[middle,0] * R)))	

#			for i in range(middle,n1):
#				for j in range(n2):
#					res[i,j] = p00*gamma/(A + a * zarg[middle - i,j])**2
			return res


#		def cs00(z):
#			zarg = z[0]
#			if(not hasattr(zarg,"__len__")):
#				if(zarg-z0[0]<0.5*(zf[0]-z0[0])):
#					return  A - a *zarg
#				return  A + a *zarg
#			n = len(zarg)
#			res = np.zeros((n,n))
#			middle = int(n/2)
#			for i in range(middle):
#				for j in range(n):
#					res[i,j] = A - a *zarg[i,j]
#			for i in range(middle,n):
#				for j in range(n):
#					res[i,j] = A + a *zarg[middle - i,j]
#			return res
	

	cs00 = lambda z: np.sqrt(gamma * p00 / rho0(z))

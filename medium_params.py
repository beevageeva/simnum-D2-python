import numpy as np
import math
from constants import gamma, z0, zf
from common import getArrayZShape
from sys import exit


p00 = 1.0
v00 = 0.0

#mediumType = "homog"
mediumType = "inhomog"  #variable density rho00 to test with wave packet

if(mediumType=="homog"):
	#rho00 = 2.0
	#rho00 = 0.5
	rho00 = 1.0
	cs00 = math.sqrt(gamma * p00 / rho00)

elif(mediumType=="inhomog"):
	#inhomogSubtype = 1
	inhomogSubtype = 2

	#rhoType = 1
	#rhoType = 2
	rhoType = 3


	if rhoType == 1 or rhoType ==2:
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
		if inhomogSubtype == 1:
			rho00 = 1.0
			rho01 = 0.01
		elif inhomogSubtype == 2:
			rho00 = 0.3
			rho01 = 1.2
			#rho01 = 1.0
			#rho00 = 1.0
			#rho01 = 2.0
		
		#rho00 = 0.3  #second exp of inhom
		#rho01 = 1.2#second exp of inhom
		ze = [0.5*(z0[0] + zf[0]),0.5*(z0[1] + zf[1])]
		#ze = [z0[0] + 0.7*(zf[0] - z0[0]),z0[1] + 0.7*(zf[1] - z0[1])]#second exp of inhom
		we = 0.4
		#we = 0.5#second exp of inhom
		#I have to apply func (argument function) before applying tanh: see initcond_soundwave
		densFunc = lambda z: 1 + np.tanh(z /we)
	
		
	
		
		def rho0(z):
			from common import getArrayZShape
			return rho00 + 0.5 * (rho01-rho00) * densFunc(densargFunc(z - getArrayZShape(ze[0], ze[1], len(z[0]))))

		cs00 = lambda z: np.sqrt(gamma * p00 / densFunc(densargFunc(z)))


	elif rhoType == 2:
		#c of form -xi + A
		rho00 = 1
		a = 3.0
		A = np.sqrt(gamma * p00 / rho00)
		def rho0(z):
			return  p00*gamma/(A + a * densargFunc(z))**2 
		cs00 =  lambda z: A + a * densargFunc(z)

	
	elif rhoType == 3:
		#c of form -xi + A
		rho00 = 1
		a = 3.0
		A = np.sqrt(gamma * p00 / rho00)
		from constants import z0,zf
		def rho0(z):
			zarg = z[0]
			if(not hasattr(zarg,"__len__")):
				if(zarg-z0[0]<0.5*(zf[0]-z0[0])):
					return p00*gamma/(A - a * zarg)**2  
				return  p00*gamma/(A + a * zarg)**2
			n = len(zarg)
			middle = int(n/2)
			res = np.zeros((n,n))
			for i in range(middle):
				for j in range(n):
					res[i,j] = p00*gamma/(A - a * zarg[i,j])**2
			for i in range(middle,n):
				for j in range(n):
					res[i,j] = p00*gamma/(A + a * zarg[i,j])**2
			return res


		def cs00(z):
			zarg = z[0]
			if(not hasattr(zarg,"__len__")):
				if(zarg-z0[0]<0.5*(zf[0]-z0[0])):
					return  A - a *zarg
				return  A + a *zarg
			n = len(zarg)
			res = np.zeros((n,n))
			middle = int(n/2)
			for i in range(middle):
				for j in range(n):
					res[i,j] = A - a *zarg[i,j]
			for i in range(middle,n):
				for j in range(n):
					res[i,j] = A + a *zarg[i,j]
			return res
	


"""
	The parameters for gaussian sound wave packet
	k0 - wavenumber
	zc - center of the gauss function
	W width of the gauss function

"""

import numpy as np
from constants import z0, zf
from math import pi,sqrt
from perturbation_params import argFunc

k0 = 60.0
#k0 = 30.0
#k0 = 50.0
#k0 = 15.0
#W = 0.05
W = 0.2
#W = 0.15

#with rhoType = 1 k0=50 W = 0.15 works
#with rhoType = 2 k0=60 W = 0.25 works ?

#W = 0.025
#W = 0.1

from medium_params import mediumType
if mediumType == "homog":
	zc = [0.5 * (z0[0]+ zf[0]), 0.5 * (z0[1] + zf[1])] #in the middle
elif mediumType == "inhomog":
	from medium_params import rhoType
	if rhoType == 4:
		zc = [z0[0] + 0.5 * (zf[0] - z0[0]), z0[1] + 0.15 * (zf[1] - z0[1])]
	else:
		zc = [z0[0] + 0.2 * (zf[0] - z0[0]), z0[1] + 0.2 * (zf[1] - z0[1])]

#zc = [0.65 * (z0[0]+ zf[0]), 0.65 * (z0[1] + zf[1])] #at the end to test reflection

def getSoundWaveGaussFunction(zc, W):
	"""
	returns the gauss function (of the envelope) - this is a gaussian wave packet!	
		Parameters
		----------
		zc : real
			center of the gauss function    
		W : real
			width of the gauss function
	returns g so that
	g(z) = -(z-zc)**2/W**2

	"""
	def gaussFunction(z):
		t2 = np.subtract(z[0],zc[0]) ** 2 + np.subtract(z[1], zc[1]) ** 2
		return np.exp(-np.divide(t2, W**2))	 
	return gaussFunction


def getSoundWaveFunction(argFunc, k0, zc, W):
	"""
	returns the resulting wave function obtained by multiplying the gauss function with a cos function:
		Parameters
		----------
		k0 : int  
			wavenumber
		zc : real
			center of the gauss function    
		W : real
			width of the gauss function
	returns f so that:
	f(z) =  g(z) / cos(2 pi k0 (z-z0)/(zf-z0) )
	with g the gauss function defined above

	"""
	def gaussPacketFunction(z):
		from common import getArrayZShape

		return np.multiply(getSoundWaveGaussFunction(zc, W)(z),  np.cos(2.0 * pi * k0 * argFunc(z-getArrayZShape(z0[0], z0[1])) ) )
	return gaussPacketFunction






def getTrajectory(z, time):
	from medium_params import cs00
	from common import getZIndex0, getZIndex1, getDz0, getDz1
	from perturbation_params import k1,k2
	if not hasattr(k1, "__len__"):
		k1 = [k1]
		k2 = [k2]
	newz = np.ones(z[0].shape)	
	tt = 0
	dt = 0.01
	dz0 = getDz0()
	dz1 = getDz1()
	ind0 = getZIndex0(zc[0])	
	ind1 = getZIndex0(zc[1])
	newz[ind0,ind1] = 0
	from medium_params import mediumType
	if mediumType == "inhomog":
		#dont use numpy gradient	
		#cs00z = cs00(z)
		#from common import derivZ0, derivZ1
		#gradCs0 = derivZ0(cs00z)
		#gradCs1 = derivZ1(cs00z)
		gradCs = np.gradient(cs00(z), getDz0(), getDz1())
	for i in range(len(k1)):
		print("k1=%e,k2=%e" % (k1[i],k2[i]))	
		lastk = [np.float128(2.0* pi *k1[i]* k0), np.float128(2.0* pi *k2[i]* k0)]
		#DO NOT SET LASTX TO zc changing lastx will change zc(same object reference)!  (do not assign:  lastx = zc)
		lastx = [zc[0],zc[1]]
		print("%e %e SET LASTX" % (lastx[0], lastx[1])	)
	
	
		stop = False	
		#np.set_printoptions(threshold='nan')
		#print("GRADIENT")
		#print(gradCs)
		lat = 2
		while(tt<time and not stop):	
			modlastk = (lastk[0]**2 + lastk[1]**2)**0.5
			indLastX0 = getZIndex0(lastx[0])	
			indLastX1 = getZIndex1(lastx[1])
			if mediumType == "homog":
				cs = cs00
				gcs0 = 0
				gcs1 = 0
			else:
				cs = cs00(lastx)
				gcs0 = gradCs[0][indLastX0,indLastX1]
				gcs1 = gradCs[1][indLastX0,indLastX1]
				#gcs0 = gradCs0[indLastX0,indLastX1]
				#gcs1 = gradCs1[indLastX0,indLastX1]
			
			#print("%e %e %d %d" % (lastx[0], lastx[1],indLastX0, indLastX1))	
	
			lastx[0]+=dt*cs*lastk[0]/modlastk
			lastx[1]+=dt*cs*lastk[1]/modlastk
			lastk[0]+=dt* np.float128(-modlastk * gcs0)
			lastk[1]+=dt*np.float128(-modlastk * gcs1)
			from os import system
			#system("echo \"%d %d %e %e %e %e %e %e  \" >> traj " % (indLastX0, indLastX1, gcs0, gcs1, lastk[0], lastk[1], dt* np.float128(-modlastk * gcs0),dt*np.float128(-modlastk * gcs1) ))
			if lastx[0] <= z0[0] + lat * dz0 or  lastx[0] >= zf[0] - lat * dz0 or  lastx[1] <= z0[1] + lat * dz0 or lastx[1] >= zf[1] - lat * dz1:
				stop = True
				print("packet out of grid traj STOP")
				
			else:
				ind0 = getZIndex0(lastx[0])	
				ind1 = getZIndex1(lastx[1])
				for i in range(ind0-lat,ind0+lat):
					for j in range(ind1-lat,ind1+lat):
						newz[i,j] = 0
			tt+=dt
	return newz
	

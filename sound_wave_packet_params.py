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

#k0 = 60.0
k0 = 15.0
W = 0.07
#W = 0.025
#W = 0.1

from medium_params import mediumType
if mediumType == "homog":
	zc = [0.5 * (z0[0]+ zf[0]), 0.5 * (z0[1] + zf[1])] #in the middle
elif mediumType == "inhomog":
	from medium_params import inhomogSubtype
	if inhomogSubtype == 1:
		zc = [z0[0] + 0.2 * (zf[0] - z0[0]), z0[1] + 0.2 * (zf[1] - z0[1])]
	
	elif inhomogSubtype == 2:
		#change 	
		k0 = 15.0 #second exp of inhom
		#zc = z0 + 0.2 * (zf - z0) homog
		#zc = z0 + 0.3 * (zf - z0)
		zc = [z0[0] + (3.0/20.0)*(zf[0] - z0[0]), z0[1] + (3.0/20.0)*(zf[1] - z0[1])  ]#second exp of inhom and first new
		#W = 0.01
		W = 0.25 #second exp of inhom

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


def getSoundWaveFunction(k0, zc, W):
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




def getSoundWaveFFTAnalytical(k0, zc, W):
	"""
		Parameters
		----------
		k0 : int  
			wavenumber
		zc : real
			center of the gauss function    
		W : real
			width of the gauss function
	
		returns analytical expression of the sound wave function defined above: [-(z-zc)**2/W**2] *  cos(2 pi k0 (z-z0)/(zf-z0) )
		calculated with mathematica 
				
	
	"""
	from constants import z0, zf
	def analyticFFT(k):
		return ((np.exp(-((pi*(k0**2*pi*W**2 - 2*k0*(k*pi*W**2 - 1j*z0 + 1j*zc)*(z0 - zf) + k*(k*pi*W**2 + (2*1j)*zc)*(z0 - zf)**2))/(z0 - zf)**2)) +  np.exp(-((pi*(k0**2*pi*W**2 + 2*k0*(k*pi*W**2 - 1j*z0 + 1j*zc)*(z0 - zf) + k*(k*pi*W**2 + (2*1j)*zc)*(z0 - zf)**2))/(z0 - zf)**2)))*sqrt(pi)*W)/2
	return analyticFFT

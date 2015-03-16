import numpy as np
import sys,math
from constants import gamma



		
			


def getGradientFunctionSym(functionType):
	try:
		if functionType == "sine":
			from sound_wave_sine_params import getSoundWaveSymDerivFunction
		elif functionType == "gauss":
			from sound_wave_gauss_params import getSoundWaveSymDerivFunction
		elif functionType == "wavepacket":
			from sound_wave_packet_params import getSoundWaveSymDerivFunction
		elif functionType == "hankel":
			from sound_wave_hankel_params import getSoundWaveSymDerivFunction
		elif functionType == "defined":
			from sound_wave_defined_params import getSoundWaveSymDerivFunction
	except ImportError:
		return None
	return getSoundWaveSymDerivFunction




def getFunctionFromType(functionType):
	try:
		if functionType == "sine":
			import sound_wave_sine_params
			from sound_wave_sine_params import phi0, k0, getSoundWaveFunction
			w = getSoundWaveFunction(k0, phi0)
		elif functionType == "gauss":
			from sound_wave_gauss_params import getSoundWaveFunction,R
			w = getSoundWaveFunction(R)
		elif functionType == "wavepacket":
			from sound_wave_packet_params import getSoundWaveFunction, k0, zc, W
			w = getSoundWaveFunction(k0, zc, W)
		elif functionType == "hankel":
			from sound_wave_hankel_params import getSoundWaveFunction
			from sound_wave_monochr_params import k
			w = getSoundWaveFunction(k)
		elif functionType == "defined":
			from sound_wave_defined_params import getSoundWaveFunction
	except ImportError:	
		print("Function type %s not defined !" % functionType)
		import sys
		sys.exit(0)
	return w

def getInitialPresRhoVel(z):
	from medium_params import p00, v00, cs00, rho00,  mediumType
	from perturbation_params import A,  functionType, waveType
	from constants import z0, zf

	w = getFunctionFromType(functionType)
	f = w(z)

	if(waveType == "lineal"):
		from perturbation_params import argFunc, velFunc


		v1 = v00 + cs00 * A* f
		vel = np.dstack(velFunc(v1))
		return {'pres': p00 + gamma * p00 * A* f  , 'rho': rho00 + rho00 *A* f , 'vel': vel }
	elif(waveType == "radial"):
		#now w is the velocity potential # gradient
		deriv = getGradientFunctionSym(functionType)
		if(deriv is None):
			#calculate numerically
			from common import getDz0, getDz1
			from math import sqrt	
			dz0 = getDz0()
			dz1 = getDz1()
			grad =  np.gradient(f, dz0, dz1 )
			fder = np.sqrt(grad[0]**2 + grad[1]**2)
		else:
			fder = deriv(np.sqrt(z[0]**2+z[1]**2))
		#print("f")	
		#print(f.shape)	
		#print("fder")	
		#print(fder.shape)	
		from math import sqrt
		v1 = (v00 + A * np.real(fder)) * sqrt(2)			
		vel = np.dstack((v1, v1))
		n = len(z[0])
		for i in range(n/2):
			for j in range(n/2):
				#change sign and assure symmetric
#				vel[n/2-i-1][n/2-j-1][0] = -vel[n/2+i][n/2+j][0]
#				vel[n/2-i-1][n/2-j-1][1] = -vel[n/2+i][n/2+j][1]
#				vel[n/2-i-1][n/2+j][0] = vel[n/2+i][n/2+j][0]
#				vel[n/2-i-1][n/2+j][1] = -vel[n/2+i][n/2+j][1]
#				vel[n/2+i][n/2-j-1][0] = -vel[n/2+i][n/2+j][0]
#				vel[n/2+i][n/2-j-1][1] = vel[n/2+i][n/2+j][1]

				#only change sign
				vel[n/2-i-1][n/2-j-1][0] *= -1
				vel[n/2-i-1][n/2-j-1][1] *= -1
				vel[n/2-i-1][n/2+j][1] *= -1
				vel[n/2+i][n/2-j-1][0] *= -1
				
		from perturbation_params import omega	
		presPert = A * rho00* omega * np.real(1j  * f)
		return {'pres': p00 + presPert , 'rho': rho00  + presPert / (cs00**2), 'vel': vel }






		
	






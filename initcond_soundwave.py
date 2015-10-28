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



def getRadialAnalitic(z, t):
	from perturbation_params import A,  functionType, waveType, argFunc, radialType
	from medium_params import p00, v00, cs00, rho00,  mediumType
	from sound_wave_monochr_params import k
	omega = cs00 * k
	radius = np.sqrt(z[0]**2 + z[1]**2)
	if radialType == "gradient":
		w = getFunctionFromType(argFunc, functionType)
		#I have already argFunc called in the function definition
		f = w(z)
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
			fder = deriv(z)
		#print("f")	
		#print(f.shape)	
		#print("fder")	
		#print(fder.shape)	
		
	
		v1 = (v00 + A * np.real(fder * np.exp(-1j * t * omega))) / radius 
		vel = np.dstack((v1 * z[0] , v1 * z[1]))
	#	n = len(z[0])
	#	for i in range(n/2):
	#		for j in range(n/2):
	#			#change sign
	#			vel[n/2+i][n/2+j][0] *=-1
	#			vel[n/2+i][n/2+j][1] *=-1
	#			#change sign and assure symmetric
	#			
	#
	#
	#			vel[n/2-i-1][n/2-j-1][0] = -vel[n/2+i][n/2+j][0]
	#			vel[n/2-i-1][n/2-j-1][1] = -vel[n/2+i][n/2+j][1]
	#			vel[n/2-i-1][n/2+j][0] = vel[n/2+i][n/2+j][0]
	#			vel[n/2-i-1][n/2+j][1] = -vel[n/2+i][n/2+j][1]
	#			vel[n/2+i][n/2-j-1][0] = -vel[n/2+i][n/2+j][0]
	#			vel[n/2+i][n/2-j-1][1] = vel[n/2+i][n/2+j][1]
	#
	#			#only change sign
	##				vel[n/2-i-1][n/2-j-1][0] *= -1
	##				vel[n/2-i-1][n/2-j-1][1] *= -1
	##				vel[n/2-i-1][n/2+j][1] *= -1
	##				vel[n/2+i][n/2-j-1][0] *= -1
				
		presPert = A * rho00* omega * np.real(1j  * f * np.exp(-1j * t * omega))
	elif radialType == "stationary":
		from math import pi
		alpha = pi / 2
		presPert = 2 * A / radius * np.cos(k*radius + alpha) * np.cos(omega * t)
		v1 = v00 -2 * A / (rho00 * omega) * (-k / radius * np.sin(k * radius + alpha) - 1/ radius**2 * np.cos(k*radius + alpha) ) * np.sin(omega * t)
		vel = np.dstack((v1 * z[0] , v1 * z[1]))
	return {'pres': p00 + presPert , 'rho': rho00  + presPert / (cs00**2), 'vel': vel }





def getFunctionFromType(argFunc, functionType):
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
			w = getSoundWaveFunction(argFunc, k0, zc, W)
		elif functionType == "wavepacket_carton":
			from sound_wave_packet_carton_params import getSoundWaveFunction, k0, zc, W
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
		print(sys.exc_info()[0])
		sys.exit(0)
	return w

def getInitialPresRhoVel(z):
	from medium_params import p00, v00, cs00, rho00,  mediumType
	from perturbation_params import A,  functionType, waveType
	from constants import z0, zf
	if(waveType == "lineal"):
		from perturbation_params import argFunc, velFunc, argType
		if argType == "2d1":
			#I have already argFunc called in the function definition
			f1 = getFunctionFromType(argFunc[0], functionType)(z) 
			f2 = getFunctionFromType(argFunc[1], functionType)(z) 
		else:
			f = getFunctionFromType(argFunc, functionType)(z)
	
		if mediumType == "homog":
			cs0 = cs00
			rho0 = rho00
		else:
			from medium_params import rho0
			cs0 = cs00(z)
			rho0 = rho0(z)
		if argType == "2d1":
			v1 = v00 + cs0 * A* f1
			v2 = v00 + cs0 * A* f2
			vfunc1 = velFunc[0](v1) 
			vfunc2 = velFunc[1](v2) 
			vel1 = np.dstack(vfunc1)
			vel2 = np.dstack(vfunc2)
			vel = vel1 + vel2
			f = f1+f2
		else:
			v1 = v00 + cs0 * A* f
			vfunc = velFunc(v1)
			vel = np.dstack(vfunc)
		return {'pres': p00 + gamma * p00 * A* f  , 'rho': rho0 + rho0 *A* f , 'vel': vel }
	elif(waveType == "radial"):
		return getRadialAnalitic(z, 0)


def getAnRhoPresVel(z, time):
	from perturbation_params import A,  functionType, waveType
	if(waveType == "radial"):
		return getRadialAnalitic(z, time)
	else:
		import sys
		print("analytical sol for other than radial not implemented")
		sys.exit(0)



		
	






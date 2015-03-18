import numpy as np
import sys,math
from constants import gamma
from medium_params import mediumType
from perturbation_params import wType



		
if wType == "pot":


	def getGradientFunction(f, z):
		try:
			from perturbation_params import symDerivW
			return symDerivW(z)
		except ImportError:
			from common import getDz0, getDz1
			dz0 = getDz0()
			dz1 = getDz1()
			return np.gradient(f,dz0, dz1 )
			


	def getAnValuesFromVelPot(z, t):
		#only for homogenous medium, otherwise exists(see sound_wave_params)
		from medium_params import p00, rho00, v00
		from perturbation_params import k, A, velFunc, w
		from common import derivZ0, derivZ1
		cs00 = math.sqrt(gamma * p00 / rho00) #only for homogenous
		f = w(z)
		omega = k * cs00
		gradW = getGradientFunction(f, z)
		np.set_printoptions(threshold='nan')
		v1 = np.real(np.multiply(gradW, np.exp(-1j * omega * t)))
		#TODO velPert = A* np.dstack(velFunc(v1, z))
		#do NOT APPLY velFunc 
		velPert = A* np.dstack((v1[0],v1[1]))
		presPert = rho00 * omega * A * np.real(1j * f * np.exp(-1j * omega * t))
		rhoPert = presPert / cs00 ** 2
		return {'pres': p00 + presPert  , 'rho': rho00 + rhoPert , 'vel': v00 + velPert } 





def getInitialPresRhoVel(z):
	from medium_params import p00, v00, cs00, rho00,  mediumType
	from perturbation_params import A,  w, wType,timesZArgW
	f =  w(z)
	if  wType == "all":
		from perturbation_params import velFunc
		from medium_params import  mediumType
		if(mediumType == "inhomog"):
			from medium_params import rho0	
			rho00 = rho0(z)
			cs00 = cs00(z)
		#initial velocity
		v1 = v00 + cs00 * A* f
		vel = np.dstack(velFunc(v1, z))
		return {'pres': p00 + gamma * p00 * A* f  , 'rho': rho00 + rho00 *A* f , 'vel': vel }
	elif wType == "pot":
		return getAnValuesFromVelPot(z, 0)



def getAnRhoPresVel(z, t):
	from perturbation_params import timesZArgW
	from medium_params import  mediumType
	if mediumType !="homog":
		print("analitical solution not implemented FOR inhomogen  medium  set plotAnalitical = False in notifier_params.py")
		sys.exit(0)
	from perturbation_params import A, wType, w, timesZArgW
	from medium_params import   p00, rho00, v00
	cs00 = math.sqrt(gamma * p00 / rho00) #only for homogenous
	#TODO chapuza
	if(timesZArgW==1):
		from perturbation_params import argType
	if(wType == "pot"):
		return getAnValuesFromVelPot(z, t)
	#no analitical solution for wave packet

	from boundary_conditions import periodicType	
	if(wType == "all" and timesZArgW == 1 and (argType == "x" or argType == "y") and periodicType == "repeat"):
		from perturbation_params import velFunc
		from common import getPeriodicXArray2
		from constants import z0_0, zf_0, z0_1, zf_1
		if(argType == "x"):
			newz0 = z[0] - cs00 * t
			z = [getPeriodicXArray2(newz0, z0_0, zf_0), z[1]]
		elif(argType == "y"):
			newz1 = z[1] - cs00 * t
			z = [z[0], getPeriodicXArray2(newz1, z0_1, zf_1)]
		f = w(z)
		v1 = v00 + cs00 * A* f
		vel = np.dstack(velFunc(v1, z))
		return {'pres': p00 + gamma * p00 * A* f  , 'rho': rho00 + rho00 *A* f , 'vel': vel }
	else:
		print("analitical solution not implemented set plotAnalitical = False in notifier_params.py")
		sys.exit(0)

		
	






import numpy as np
import math
from constants import gamma, z0_0, zf_0, z0_1, zf_1
from common import getArrayZShape
from sys import exit


p00 = 1.0
v00 = 0.0

mediumType = "homog"
#mediumType = "inhomog"  #variable density rho00 to test with wave packet

if(mediumType=="homog"):
	rho00 = 1.0
	cs00 = math.sqrt(gamma * p00 / rho00)

elif(mediumType=="inhomog"):
	rho00 = 1.0
	rho01 = 0.01
	#rho00 = 0.3  #second exp of inhom
	#rho01 = 1.2#second exp of inhom
	ze = [0.5*(z0_0 + zf_0),0.5*(z0_1 + zf_1)]
	#ze = [z0_0 + 0.7*(zf_0 - z0_0),z0_1 + 0.7*(zf_1 - z0_1)]#second exp of inhom
	we = 0.4
	#we = 0.5#second exp of inhom
	#densFunc = lambda z: 1 + np.tanh((z - getArrayZShape(ze[0], ze[1], len(z[0])))/we)
	#I have to apply func (argument function) before applying tanh: see initcond_soundwave
	densFunc = lambda z: 1 + np.tanh(z /we)
	from common import getArrayZShape
	from perturbation_params import densargfunc	
	def rho00(z):
		return rho00 + 0.5 * (rho01-rho00) * densFunc(densargfunc(z - getArrayZShape(ze[0], ze[1], len(z[0]))))
	
	cs00 = lambda(z): np.sqrt(gamma * p00 / densFunc(z))


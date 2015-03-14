import numpy as np
from math import pi
from constants import z0, zf
from perturbation_params import argFunc

phi0 = pi / 6.0 
k0 = 1



def getSoundWaveFunction(k0, phi0):
		def sinFunction(z):
			from common import getArrayZShape
			#return np.sin(2.0 * pi * k0 * argFunc(z - getArrayZShape(z0[0], z0[1]))/wl + phi0 )
			return np.sin(2.0 * pi * k0 *  argFunc(z - getArrayZShape(z0[0], z0[1])) + phi0 )
			#return np.sin(2.0 * pi * k0 * argFunc(z) + phi0 )
			#np.sin(np.multiply((2.0 * math.pi),argFunc(z-z0)) + phi0 )
			#return np.sin(np.multiply((2.0 * pi),argFunc(z)) + phi0 )
		return sinFunction

import numpy as np
from math import pi
from constants import z0, zf
from perturbation_params import argFunc




phi0 = pi / 6.0 
k0 = 1

"""
	returns the sondWaveFunction
	Note I have to apply here argFunc from perturbation_params becuase soemtimes I have to apply it to some other arg as z0
"""

def getSoundWaveFunction(k0, phi0):
		def sinFunction(z):
			from common import getArrayZShape
			return np.sin(2.0 * pi * k0 *  argFunc(z - getArrayZShape(z0[0], z0[1])) + phi0 )
		return sinFunction

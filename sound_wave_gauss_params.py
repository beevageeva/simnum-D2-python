import numpy as np

from constants import z0, zf
from perturbation_params import argFunc


R = 0.05
def getSoundWaveFunction(R):
	def w(z):
		from common import getArrayZShape
		middlex = 0.5 * (z0[0] + zf[0])	
		middley = 0.5 * (z0[1] + zf[1])	
		return np.exp(-argFunc(z - getArrayZShape(middlex, middley))**2 / R)
	return w;

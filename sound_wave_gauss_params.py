import numpy as np

from constants import z0, zf

R = 0.05
def getSoundWaveFunction(R):
	def w(z):
		return np.exp(-(z - 0.5 *(z0 + zf))**2 / R)
	return w;

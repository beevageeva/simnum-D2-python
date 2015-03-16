from scipy.special import hankel1
import numpy as np

from sound_wave_monochr_params import k


smoothInterp = False
#smoothInterp = True
withWindow = False
#withWindow = True
useAnDer = True
#useAnDer = False

"""
	returns the sondWaveFunction
	Note I have to apply here argFunc from perturbation_params becuase soemtimes I have to apply it to some other arg as z0
"""
def getSoundWaveFunction(k):
	from perturbation_params import argFunc
	def resFunc(z): 
		res =  hankel1(0, k * argFunc(z))
		if smoothInterp:
			smoothInterp(z, res)
		if withWindow:
			#multiply by a window function
			windowFunction1D = np.blackman(len(z[0]))
			windowFunction2D = np.outer(windowFunction1D, windowFunction1D)
			return res * windowFunction2D
		else:
			return res
	return resFunc

#if smooth interp or multiplied by a window function this will not be the gradient anymore
#if not smoothInterp and not withWindow:
if useAnDer:
	"""
		returns the nalytic derivative of soundWaveFunction
		Note I have to apply here argFunc from perturbation_params
	"""
	def getSoundWaveSymDerivFunction(z):
		from perturbation_params import argFunc
		res =  hankel1(1, k * argFunc(z))
		if smoothInterp:
			smoothInterp(z, res)
		if withWindow:
			#multiply by a window function
			windowFunction1D = np.blackman(len(z[0]))
			windowFunction2D = np.outer(windowFunction1D, windowFunction1D)
			res *= windowFunction2D
		return -k * res



def smoothInterp(z, res):
	from scipy.interpolate import CloughTocher2DInterpolator
	RI = 10 #mask a circle around 0,0 -> assumes a centered domain with radius RI and sets  the values here by intterpolation
	#n = len(z[0]) - 1
	#y,x=np.ogrid[-n / 2: n/2 + 1, -n / 2: n/2 + 1]
	n = len(z[0]) 
	y,x=np.ogrid[-n / 2: n/2 , -n / 2: n/2 ]
	mask = x**2+y**2 < RI**2
	nMask = ~mask
	acl = np.dstack((z[0][nMask], z[1][nMask]))
	acl = acl[0] #TODO why
	#real and imag comp separately: some interpolation functions do not handle complex values but CloughTocher2DInterpolator does
#			freal = CloughTocher2DInterpolator(acl, np.real(res[nMask]))
#			resIntReal = freal(z[0][mask], z[1][mask])
#			fimag = CloughTocher2DInterpolator(acl, np.imag(res[nMask]))
#			resIntImag = fimag(z[0][mask], z[1][mask])
#			res[mask] = resIntReal + 1j * resIntImag
	f = CloughTocher2DInterpolator(acl, res[nMask])
	#print("n=%d" % n)
	res[mask] = f(z[0][mask], z[1][mask])
	for i in range(n/2):
		for j in range(n/2):
			res[n/2-i-1][n/2+j] = res[n/2+i][n/2+j]
			res[n/2-i-1][n/2-j-1] = res[n/2+i][n/2+j]
			res[n/2+i][n/2-j-1] = res[n/2+i][n/2+j]

def w1(z, z0):
	#r = k * func(z-z0)
	r = k * func(z)
	res = hankel1(0, r)
	smoothInterp(z, res)
	#multiply by a window function
	#n = len(z[0]) - 1
	#windowFunction1D = np.blackman(n+1)
	#windowFunction2D = np.outer(windowFunction1D, windowFunction1D)
	#return res * windowFunction2D
	return res

			


		#uncomment the following in order to user analitic deriv of w, otherwise it will be numerically calculated using numpy.gradient: see initcont_soundwave
		#take care when function is multiplied by a window function the following is NOT the gradient
#		def symDerivW(z):
#			d1Arg = -k * hankel1(1, k * func(z))
#			smoothInterp(z, d1Arg)
#			#multiply by a window function
#			n = len(z[0]) - 1
#			windowFunction1D = np.blackman(n+1)
#			windowFunction2D = np.outer(windowFunction1D, windowFunction1D)
#			d1Arg *= windowFunction2D
#			return velFunc([d1Arg, d1Arg], z)

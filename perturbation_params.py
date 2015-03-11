import numpy as np
import math
from constants import gamma, z0_0, zf_0, z0_1, zf_1
from common import getArrayZShape
from sys import exit


from medium_params import mediumType

A = 3.0 * 10.0 ** (-4)
#A = 5.0 * 10.0 ** (-2)


wType = "all"
#wType = "pot"

if(wType=="pot" and mediumType == "inhomog"):
	print("INVALID config: wType=pot and inhomog medium")
	sys.exit(0)


timesZArgW = 1 #1(sine, gauss, hankel test - with wType = "pot") or 2(wave packet)
#timesZArgW = 2 

if(timesZArgW == 1):
	functionType = "sine" 
	#functionType = "gauss" 
	#functionType = "hankel" 
	argType = "x"
	#argType = "y"
	#argType = "r"
	#argType = "d1" 
	

	if(argType == "x"):
		nx = 1.0
		wl = zf_0 - z0_0
		k1 = nx / (zf_0 - z0_0)
		velFunc1 = lambda x1, x2, z: (x1,np.zeros(x2.shape))
		func = lambda z: k1 * z[0] 
	elif argType == "y":
		wl = zf_1 - z0_1
		ny = 1
		k2 = ny / (zf_1 - z0_1)
		velFunc1 = lambda x1, x2, z: (np.zeros(x1.shape),x2)
		func = lambda z: k2 * z[1] 
	elif argType == "d1":
		wl =  math.sqrt((zf_0 - z0_0)**2 + (zf_1 - z0_1)**2)
		nx= 2.0
		ny = 3.0	
		k1 = nx/(zf_0 - z0_0)
		k2 = ny/(zf_1 - z0_1)
		modk = math.sqrt(k1**2 + k2**2)
		velFunc1 = lambda x1 ,x2, z: (k1/modk * x1 ,k2/modk * x2)
		func = lambda z: z[0] * k1 + z[1] * k2 
	elif argType == "r":
		wl =  math.sqrt((zf_0 - z0_0)**2 + (zf_1 - z0_1)**2)
		def velFunc1(x1, x2, z):
			r = np.sqrt(z[0]**2 + z[1]**2)
			return x1* (z[0]/r) , x2 * (z[1] / r) 
		func = lambda z: np.sqrt(z[0] **2 +  z[1] **2)

	if wType=="pot":
		k0 = 3
		k = k0 * 2 * math.pi / wl
		velFunc = lambda x, z: velFunc1(x[0], x[1], z)
	else:	
		velFunc = lambda x, z: velFunc1(x, x, z)
	
	
	if functionType == "sine":
		phi0 = math.pi / 6.0 
		#sine function
		def w1(z,z0):
			return np.sin(np.multiply((2.0 * math.pi),func(z-z0)) + phi0 )
	
	elif functionType == "gauss":
		R = 0.05
		def w1(z,z0):
			if(argType == "r"):
				b = 0
			else:
				b = 0.5 * wl
			return np.exp(-(func(z -z0) - b)**2 / R) 

	elif functionType == "hankel":
		if(argType!="r" or wType!="pot"):
			print("argType=%s must be r for hankel and wType=%s must be pot" % (argType, wType))
			exit(0)
		from scipy.special import hankel1



		def smoothInterp(z, res):
			from scipy.interpolate import CloughTocher2DInterpolator
			RI = 1 #mask a circle around 0,0 -> assumes a centered domain with radius RI and sets  the values here by intterpolation
			n = len(z[0]) - 1
			y,x=np.ogrid[-n / 2: n/2 + 1, -n / 2: n/2 + 1]
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
			res[mask] = f(z[0][mask], z[1][mask])


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
		
	def w(z):
		n = len(z[0])
		if(argType=="r"):
			z0 = np.zeros((n,n))
		else:
			z0 = getArrayZShape(z0_0,z0_1,n)
		return w1(z, z0)

	#TODO: inhomegen medium in this case??
	if mediumType == "inhomog":
		print("!!!!mediumType = INHOMOG and timesZArgW=1(not in wave packet case)!!!")
		densargfunc = func

#wave packet
elif(timesZArgW == 2):
	kf = [2.0 * math.pi/ (zf_0 - z0_0), 2.0 * math.pi/ (zf_1 - z0_1)]
	k0 = [60.0,60.0]
	#k0 = [15.0,15.0] #second exp of inhom
	k = [k0[0] * kf[0], k0[1] * kf[1]]
	zc = [0.5 * (z0_0+ zf_0), 0.5 * (z0_1 + zf_1)] #in the middle
	#zc = [z0_0 + (3.0/20.0)*(zf_0 - z0_0),z0_1 + (3.0/20.0)*(zf_1 - z0_1) ] #second exp of inhom
	W = 0.05
	#W = 0.25#second exp of inhom
	def w(z, nwav=k0):
		n = len(z[0])
		#z0 = getArrayZShape(z0_0,z0_1,n)
		#direct operations
		#k0Reshaped = getArrayZShape(k0[0],k0[1],n)
		#zcReshaped = getArrayZShape(zc[0],zc[1],n)
		#t1 = np.subtract(z, zcReshaped)
		#t2 = t1[0,:,:] ** 2 + t1[1,:,:] ** 2
		#direct operations:
		t2 = np.subtract(z[0],zc[0]) ** 2 + np.subtract(z[1], zc[1]) ** 2
		#direct operations:
		#t3 = np.subtract(z, z0)
		#return np.multiply(np.exp(-np.divide(t2, W**2)), np.cos(t3[0] * k0[0] + t3[1] * k0[1]))
		#direct operations:
		sumt1 = np.multiply(np.subtract(z[0], z0_0), k[0])
		sumt2 = np.multiply(np.subtract(z[1], z0_1), k[1])
		return np.multiply(np.exp(-np.divide(t2, W**2)), np.cos(np.add(sumt1, sumt2)))

	modk = math.sqrt(k0[0] ** 2 + k0[1]**2)
	velFunc = lambda x, z: ((k0[0]/modk) * x , (k0[1]/modk) * x)
	if mediumType == "inhomog":
		densargfunc = lambda z: k0[0]/modk  * z[0] + k0[1] /modk * z[1] #constant in a dir perp to k
		#densargfunc = lambda z: z[0] #horrizontal
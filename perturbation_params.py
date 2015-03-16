import numpy as np
import math
from constants import gamma, z0, zf
from common import getArrayZShape
from sys import exit


from medium_params import mediumType

A = 3.0 * 10.0 ** (-4)
#A = 5.0 * 10.0 ** (-2)

waveType = "lineal"
#waveType = "radial"
functionType = "wavepacket"
#functionType = "sine" 
#functionType = "gauss" 
#functionType = "hankel" 




if waveType == "lineal":	
	#argType = "x"
	#argType = "y"
	argType = "d1" 
	

	if(argType == "x"):
		nx = 1.0
		from common import getDz0
		wl = zf[0] - z0[0] + getDz0()
		k1 = nx / wl
		#k1 = nx 
		argFunc = lambda x: k1 * x[0]
		velFunc = lambda x: (x,np.zeros(x.shape))
		
	elif argType == "y":
		wl = zf[1] - z0[1]
		ny = 3.0
		k2 = ny / wl
		#k2 = ny
		argFunc = lambda x: k2 * x[1]
		velFunc = lambda x: (np.zeros(x.shape),x)
	elif argType == "d1":
		wl =  math.sqrt((zf[0] - z0[0])**2 + (zf[1] - z0[1])**2)
		nx= 2.0
		ny = 3.0	
		k1 = nx/(zf[0] - z0[0])
		k2 = ny/(zf[1] - z0[1])
		#k1 = nx
		#k2 = ny
		modk = math.sqrt(k1**2 + k2**2)
		argFunc = lambda x: k1 * x[0] + k2 * x[1]
		velFunc = lambda x: (k1/modk * x ,k2/modk * x)


elif waveType == "radial":
	argFunc = lambda x: np.sqrt(x[0]**2+x[1]**2)
	wl = zf[1] - z0[1]
	from sound_wave_monochr_params import k
	if mediumType == "inhomog":
		print("radial not impl for inhomog")
		exit(0)
	from medium_params import cs00
	omega = cs00 * k




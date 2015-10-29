notifierType = "visual"
#nstepsPlot = 100
nstepsPlot = 1

from constants import z0, zf

#second element in array is whether to calculate maxpoint velocity: I don't know what happens when multiple points for the max
#TODO max speed does not show expected value
#plots = {"3d": False, "dim0": [0.5 * (zf[0] + z0[0]), False], "dim1": [0.5 * (zf[1] + z0[1]), False], "line" : [lambda x,y : 2*x-3*y==0, [0.5 * (zf[0] + z0[0]),[0.5 * (zf[1] + z0[1])] False],  "color": True}

#use the following with the wave packet (no test for functionType)
from sound_wave_packet_params import zc
plots = {"3d": False, "dim0": [zc[0], False], "dim1": [zc[1], False],  "color": True}
from perturbation_params import waveType
if waveType == "lineal":
	from perturbation_params import argType
	if argType == "d1":
		from perturbation_params import nx, ny
		plots["line"] = [lambda x,y : nx*x-ny*y==0, zc, False] 	

#plots = {"3d": True, "dim0": [0.5 * (zf[0] + z0[0]), False], "dim1": [0.5 * (zf[1] + z0[1]), False], "color": True}
#plots = {"3d": False, "color": True}


#use it in wave packet
plotPresFFT = False
#plotPresFFT = True
#plotVelFFT = True
plotVelFFT = False

plotAnalitical = True
#plotAnalitical = False

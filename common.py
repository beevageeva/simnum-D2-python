from constants import nint, z0, zf
import numpy as np


def getDz0():
	return float(zf[0] - z0[0]) / nint

def getDz1():
	return float(zf[1] - z0[1]) / nint

def getZArray():
	dz0 = getDz0()
	dz1 = getDz1()
	a = np.linspace(z0[0]-0.5 *dz0, zf[0]+0.5*dz0, nint+2)
	b = np.linspace(z0[1]-0.5 *dz1, zf[1]+0.5*dz1, nint+2)
	return np.meshgrid(a, b, indexing='ij')


def getArrayZShape(x0,y0,n=nint+2):
	return np.array([[[x0,]*n,]*n,[[y0,]*n,]*n] )
	

def getPeriodicX(xval, a, b):
	p = float(b - a)
	k = int((xval-a)/p)
	res = xval - k * p
	if(res < a):
		res+=p
	if(res > b):
		res-=p
	return res


def getPeriodicXArray(xarray, a, b):
	res = []	
	for xval in xarray:
		res.append(getPeriodicX(xval, a, b))
	return np.array(res)

#xarray has 2 dim
def getPeriodicXArray2(xarray, a, b):
	res = np.zeros(xarray.shape)
	for i in range(0, xarray.shape[0]):
		for j in range(0, xarray.shape[1]):
			res[i][j] = getPeriodicX(xarray[i][j], a, b)
	return res

def getPeriodicX2(xval, a, b):
	p = b - a
	while (xval < a):
		xval+=p
	while (xval > b):
		xval-=p
	return xval

def getZIndex0(z):
	return int( float(nint+1)*(float(z) - z0[0])/(zf[0] - z0[0]) )

def getZIndex1(z):
	return int( float(nint+1)*(float(z) - z0[1])/(zf[1] - z0[1]) )


def getSpeedPeriodic(newVal, oldVal, z0, zf, dt):
	if newVal < oldVal:
		newVal +=  zf - z0	
	return (newVal - oldVal) / dt

def getSpeedPeriodic0(newVal, oldVal, dt):
	return getSpeedPeriodic(newVal, oldVal, z0[0], zf[0], dt)
	
def getSpeedPeriodic1(newVal, oldVal, dt):
	return getSpeedPeriodic(newVal, oldVal, z0[1], zf[1], dt)


#creates an output directory called out_0, out_1, ... the first that does not exists
def createFolder(dirname_base="out"):
	import os
	dirExists = True
	i = 0
	dirname = "%s_%i" % (dirname_base, i)
	while os.path.exists(dirname):
		i +=1
		dirname = "%s_%i" % (dirname_base, i)
	os.mkdir(dirname)
	return dirname

#centered
def derivZ0(f):
	dx = getDz0()
	n = len(f)
	res = np.zeros((n, n), dtype=complex )
	for i in range(0,n):	
		for j in range(1,n-1):	
			#centered
			res[i,j] = complex(f[i+1,j] - f[i-1,j])/complex(2 * dx)
		res[0,i] = res[1,i]
		res[n-1,i] = res[n-2,i]
	return res	

def derivZ1(f):
	dx = getDz1()
	n = len(f)
	res = np.zeros((n, n), dtype=complex )
	for j in range(0,n):	
		for i in range(1,n-1):	
			#centered
			res[i][j] = complex(f[i,j+1] - f[i,j-1])/complex(2 * dx)
		res[j,0] = res[j,1]
		res[j,n-1] = res[j,n-2]
	return res	




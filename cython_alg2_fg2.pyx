import numpy as np
cimport numpy as np

import cython
cimport cython

DTYPE = np.float
ctypedef np.float_t DTYPE_t




@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef calc_interm_u_array_2d(np.ndarray[DTYPE_t, ndim=2] resx, np.ndarray[DTYPE_t, ndim=2] resy, np.ndarray[DTYPE_t, ndim=2] u, np.ndarray[DTYPE_t, ndim=3] f, int nint, float dz0, float dz1, float dt):
	cdef int i,j
	for i in range(0, nint+1):
		for j in range(0, nint+1):
			resx[i,j] = 0.5 * (u[i+1,j] + u[i,j]) - dt * (0.25*(f[i,j+1,1] - f[i,j,1] + f[i+1,j+1,1] - f[i+1,j,1])/dz1 +0.5*(f[i+1,j,0] - f[i,j,0]) / dz0)
			resy[i,j] = 0.5 * (u[i,j] + u[i,j+1]) - dt * (0.25*(f[i+1,j,0] - f[i,j,0] + f[i+1,j+1,0] - f[i,j+1,0])/dz0 +0.5*(f[i,j+1,1] - f[i,j,1]) / dz1)
	for j in range(0, nint+1):
		resx[j,nint+1] = 0.5 * (u[j,nint+1] + u[j+1,nint+1]) - dt * (0.25*(f[j, nint+1,1] - f[j,nint,1] + f[j+1,nint+1,1] - f[j+1,nint,1])/dz1 +0.5*(f[j+1,nint+1,0] - f[j,nint+1,0]) / dz0)
		resy[nint+1,j] = 0.5 * (u[nint+1,j] + u[nint+1,j+1]) - dt * (0.25*(f[nint+1,j,0] - f[nint+1,j,0] + f[nint+1,j+1,0] - f[nint,j+1,0])/dz0 +0.5*(f[nint+1,j+1,1] - f[nint+1,j,1]) / dz1)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef calc_interm_u_array_3d(np.ndarray[DTYPE_t, ndim=3] resx, np.ndarray[DTYPE_t, ndim=3] resy, np.ndarray[DTYPE_t, ndim=3] u, np.ndarray[DTYPE_t, ndim=3] f, int nint, float dz0, float dz1, float dt):
	cdef int i,j
	for i in range(0, nint+1):
		for j in range(0, nint+1):
			resx[i,j,0] = 0.5 * (u[i+1,j,0] + u[i,j,0]) - dt * (0.25*(f[i,j+1,1] - f[i,j,1] + f[i+1,j+1,1] - f[i+1,j,1])/dz1 +0.5*(f[i+1,j,0] - f[i,j,0]) / dz0)
			resx[i,j,1] = 0.5 * (u[i+1,j,1] + u[i,j,1]) - dt * (0.25*(f[i,j+1,2] - f[i,j,2] + f[i+1,j+1,2] - f[i+1,j,2])/dz1 +0.5*(f[i+1,j,1] - f[i,j,1]) / dz0)
			resy[i,j,0] = 0.5 * (u[i,j,0] + u[i,j+1,0]) - dt * (0.25*(f[i+1,j,0] - f[i,j,0] + f[i+1,j+1,0] - f[i,j+1,0])/dz0 +0.5*(f[i,j+1,1] - f[i,j,1]) / dz1)
			resy[i,j,1] = 0.5 * (u[i,j,1] + u[i,j+1,1]) - dt * (0.25*(f[i+1,j,1] - f[i,j,1] + f[i+1,j+1,1] - f[i,j+1,1])/dz0 +0.5*(f[i,j+1,2] - f[i,j,2]) / dz1)
	for j in range(0, nint+1):
		resx[j,nint+1,0] = 0.5 * (u[j,nint+1,0] + u[j+1,nint+1,0]) - dt * (0.25*(f[j, nint+1,1] - f[j,nint,1] + f[j+1,nint+1,1] - f[j+1,nint,1])/dz1 +0.5*(f[j+1,nint+1,0] - f[j,nint+1,0]) / dz0)
		resx[j,nint+1,1] = 0.5 * (u[j,nint+1,1] + u[j+1,nint+1,1]) - dt * (0.25*(f[j, nint+1,2] - f[j,nint,2] + f[j+1,nint+1,2] - f[j+1,nint,2])/dz1 +0.5*(f[j+1,nint+1,1] - f[j,nint+1,1]) / dz0)
		resy[nint+1,j,0] = 0.5 * (u[nint+1,j,0] + u[nint+1,j+1,0]) - dt * (0.25*(f[nint+1,j,0] - f[nint+1,j,0] + f[nint+1,j+1,0] - f[nint,j+1,0])/dz0 +0.5*(f[nint+1,j+1,1] - f[nint+1,j,1]) / dz1)
		resy[nint+1,j,1] = 0.5 * (u[nint+1,j,1] + u[nint+1,j+1,1]) - dt * (0.25*(f[nint+1,j,1] - f[nint+1,j,1] + f[nint+1,j+1,1] - f[nint,j+1,1])/dz0 +0.5*(f[nint+1,j+1,2] - f[nint+1,j,2]) / dz1)



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef calc_final_u_array_2d(np.ndarray[DTYPE_t, ndim=2] res, np.ndarray[DTYPE_t, ndim=2] u, np.ndarray[DTYPE_t, ndim=3] intermF0, np.ndarray[DTYPE_t, ndim=3] intermF1,int n, float dz0, float dz1, float dt, int skip):
	cdef int i,j
	for i in range(0, n):
		for j in range(0, n):
			res[i,j] = u[i+skip,j+skip] - dt * ((intermF1[i+1,j,0] - intermF1[i,j,0])/dz0 + (intermF0[i,j+1,1] - intermF0[i,j,1])/dz0)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef calc_final_u_array_3d(np.ndarray[DTYPE_t, ndim=3] res, np.ndarray[DTYPE_t, ndim=3] u, np.ndarray[DTYPE_t, ndim=3] intermF0, np.ndarray[DTYPE_t, ndim=3] intermF1,int n, float dz0, float dz1, float dt, int skip):
	cdef int i,j
	for i in range(0, n):
		for j in range(0, n):
			res[i,j,0] = u[i+skip,j+skip,0] - dt * ((intermF1[i+1,j,0] - intermF1[i,j,0])/dz0 + (intermF0[i,j+1,1] - intermF0[i,j,1])/dz0)
			res[i,j,1] = u[i+skip,j+skip,0] - dt * ((intermF1[i+1,j,1] - intermF1[i,j,1])/dz0 + (intermF0[i,j+1,2] - intermF0[i,j,2])/dz0)


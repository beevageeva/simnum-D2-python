import numpy as np
cimport numpy as np

import cython
cimport cython

DTYPE = np.float
ctypedef np.float_t DTYPE_t




@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef calc_interm_u_array_2d(np.ndarray[DTYPE_t, ndim=2] res,  np.ndarray[DTYPE_t, ndim=2] u, np.ndarray[DTYPE_t, ndim=3] f, int nint1, int nint2, float dz0, float dz1, float dt):
	cdef int i,j
	for i in range(1, nint1+2):
		for j in range(1, nint2+2):
			#points displaced right +1 
			res[i-1,j-1]  = 0.25 * (u[i-1,j-1] + u[i-1,j] + u[i,j-1] + u[i,j]) - 0.25 * dt  * ((f[i,j,0] - f[i-1,j,0] + f[i,j-1,0] - f[i-1,j-1,0])/dz0 + (f[i,j,1] - f[i,j-1,1]+f[i-1,j,1] - f[i-1,j-1,1]) / dz1)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef calc_interm_u_array_3d(np.ndarray[DTYPE_t, ndim=3] res,  np.ndarray[DTYPE_t, ndim=3] u, np.ndarray[DTYPE_t, ndim=3] f, int nint1, int nint2, float dz0, float dz1, float dt):
	cdef int i,j
	for i in range(1, nint1+2):
		for j in range(1, nint2+2):
		#points displaced right +1 
			res[i-1,j-1,0] = 0.25 * (u[i-1,j-1,0] + u[i-1,j,0] + u[i,j-1,0] + u[i,j,0]) - 0.25 * dt  * ((f[i,j,0] - f[i-1,j,0]+f[i,j-1,0] - f[i-1,j-1,0]) / dz0 + (f[i,j,1] - f[i,j-1,1] + f[i-1,j,1] - f[i-1,j-1,1])/dz1 )
			res[i-1,j-1,1]  = 0.25 * (u[i-1,j-1,1] + u[i-1,j,1] + u[i,j-1,1] + u[i,j,1]) - 0.25 * dt * ((f[i,j,1] - f[i-1,j,1]+f[i,j-1,1] - f[i-1,j-1,1]) / dz0 + (f[i,j,2] - f[i,j-1,2] + f[i-1,j,2] - f[i-1,j-1,2])/dz1)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef calc_final_u_array_2d(np.ndarray[DTYPE_t, ndim=2] res, np.ndarray[DTYPE_t, ndim=2] u, np.ndarray[DTYPE_t, ndim=3] intermF, int n1, int n2, float dz0, float dz1, float dt, int skip):
	cdef int i,j
	for i in range(0, n1):
		for j in range(0, n2):
			res[i,j] = u[i+skip,j+skip] - dt  * 0.5 * ((intermF[i+1,j,0] - intermF[i,j,0] + intermF[i+1,j+1,0] - intermF[i,j+1,0])/dz0 + (intermF[i,j+1,1] - intermF[i,j,1] + intermF[i+1,j+1,1] - intermF[i+1,j,1])/dz1)

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
cpdef calc_final_u_array_3d(np.ndarray[DTYPE_t, ndim=3] res, np.ndarray[DTYPE_t, ndim=3] u, np.ndarray[DTYPE_t, ndim=3] intermF, int n1, int n2, float dz0, float dz1, float dt, int skip):
	cdef int i,j
	for i in range(0, n1):
		for j in range(0, n2):
			res[i,j,0]  = u[i+skip,j+skip,0] - dt * 0.5 * ((intermF[i+1,j,0] - intermF[i,j,0] + intermF[i+1,j+1,0] - intermF[i,j+1,0])/dz0 + (intermF[i,j+1,1] - intermF[i,j,1] + intermF[i+1,j+1,1] - intermF[i+1,j,1]) / dz1)
			res[i,j,1]  = u[i+skip,j+skip,1] - dt  * 0.5* ((intermF[i+1,j,1] - intermF[i,j,1] + intermF[i+1,j+1,1] - intermF[i,j+1,1])/dz0 + (intermF[i,j+1,2] - intermF[i,j,2] + intermF[i+1,j+1,2] - intermF[i+1,j,2])/dz1)

